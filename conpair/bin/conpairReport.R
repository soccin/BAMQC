#!/usr/bin/env Rscript
##
## Report for cohort-wide Conpair. Reads the tables conpairCohort.py wrote
## into <OutDir>, groups samples by genotype, flags identity and
## contamination problems, and writes summary tables, a text report and a
## PDF.
##
## Run with --help for usage. Exits 0 whenever the inputs could be read:
## the report is advisory and leaves the decision to a person.
##

suppressPackageStartupMessages({
    library(tidyverse)
    library(fs)
    library(glue)
    library(ComplexHeatmap)
    library(patchwork)
})

##
## Tool and repo directories, from this script's own location in bin/.
##
scriptPath <- commandArgs(trailingOnly = FALSE) |>
    str_subset("^--file=") |>
    str_remove("^--file=")
toolDir <- path_dir(path_dir(path_real(scriptPath)))
repoRoot <- path_dir(toolDir)
source(path(repoRoot, "R", "palettes.R"))
VERSION <- read_lines(path(toolDir, "VERSION"))[1]

##
## Arguments.
##
args <- commandArgs(trailingOnly = TRUE)

getOpt <- function(flag, fallback) {
    i <- which(args == flag)
    if (length(i) > 0 && length(args) > i[1]) args[i[1] + 1] else fallback
}

usage <- function() {
    glue("
conpairReport.R {VERSION} -- report for cohort-wide Conpair

Usage:
  Rscript conpairReport.R <OutDir> [options]

<OutDir> is the output directory of runConpairCohort.sh. Reads
manifest.tsv and the conpair_*.csv tables written by conpairCohort.py.

Options (all in percent except --min-markers):
  --match <x>        Identity at or above this means the same individual.
                     Default: the value conpairCohort.py used, from
                     conpair_run.csv, so partners and groups agree.
  --related <x>      Identity at or above this between samples labeled as
                     different individuals is a PARTIAL_MATCH warning.
                     Default: 80
  --min-markers <n>  A concordance based on fewer markers is not used.
                     Default: the value conpairCohort.py used.
  --contam-warn <x>  Contamination warning threshold. Default: 2
  --contam-fail <x>  Contamination failure threshold. Default: 5
  -h, --help         Show this message and exit.

Writes into <OutDir>:
  conpair_summary.csv   one row per sample, worst QC first
  conpair_pairs.csv     one row per sample pair, flagged pairs first
  conpair_tables.xlsx   every conpair_*.csv table as a sheet
  conpair_report.txt   text report, also printed to stdout
  conpair_report.pdf    overview, flagged samples, concordance matrix, one
                        small matrix per individual, contamination per
                        sample, glossary
") |> cat("\n")
    quit(status = 0)
}

if (length(args) == 0 || any(args %in% c("-h", "--help")) || str_starts(args[1], "-")) {
    usage()
}
outDir <- args[1]

readTable <- function(name) {
    read_csv(path(outDir, name), na = "NA", show_col_types = FALSE)
}

runTable <- read_csv(path(outDir, "conpair_run.csv"),
                     col_types = cols(.default = col_character()),
                     show_col_types = FALSE)
run <- deframe(runTable)

matchPct <- as.numeric(getOpt("--match", run[["match"]]))
minMarkers <- as.numeric(getOpt("--min-markers", run[["minMarkers"]]))
relatedPct <- as.numeric(getOpt("--related", "80"))
contamWarn <- as.numeric(getOpt("--contam-warn", "2"))
contamFail <- as.numeric(getOpt("--contam-fail", "5"))

##
## Inputs. The manifest is read with only "" as missing, so a sample can
## never turn into NA.
##
manifest <- read_tsv(path(outDir, "manifest.tsv"),
                     col_types = cols(.default = col_character()),
                     na = "", show_col_types = FALSE)
sampleOrder <- manifest$sample
hasIndividual <- any(!is.na(manifest$individual))

samples <- readTable("conpair_samples.csv")
contam <- readTable("conpair_contamination.csv") |>
    mutate(partner = as.character(partner))

##
## Directed concordance. A value resting on fewer than minMarkers markers
## is dropped rather than trusted.
##
concordance <- readTable("conpair_concordance.csv")
directed <- concordance |>
    transmute(
        reference, query,
        concHom = if_else(nMarkersHom >= minMarkers, concordanceHom, NA_real_),
        nMarkersHom,
        concAll = if_else(nMarkersAll >= minMarkers, concordanceAll, NA_real_),
        nMarkersAll
    )

##
## Pairs. Identity is the larger of the two directed concordances: LOH in
## one sample lowers only the direction where that sample is the
## reference, so the larger value is the one LOH has not affected.
##
pairs <- directed |>
    filter(match(reference, sampleOrder) < match(query, sampleOrder)) |>
    select(sampleA = reference, sampleB = query,
           concordanceAB = concHom, nMarkersAB = nMarkersHom,
           concordanceAll = concAll, nMarkersAll) |>
    left_join(
        directed |>
            select(sampleA = query, sampleB = reference,
                   concordanceBA = concHom, nMarkersBA = nMarkersHom),
        by = c("sampleA", "sampleB")
    ) |>
    mutate(
        identity = pmax(concordanceAB, concordanceBA, na.rm = TRUE),
        asymmetry = round(abs(concordanceAB - concordanceBA), 2)
    )

##
## Genotype groups: connected components of the graph whose edges are
## pairs at or above the match threshold. Numbered in manifest order.
##
graph <- igraph::graph_from_data_frame(
    pairs |> filter(identity >= matchPct) |> select(sampleA, sampleB),
    directed = FALSE,
    vertices = tibble(name = sampleOrder)
)
membership <- igraph::components(graph)$membership
groups <- tibble(sample = names(membership), component = membership) |>
    arrange(match(sample, sampleOrder)) |>
    mutate(groupIndex = match(component, unique(component))) |>
    mutate(genotypeGroup = str_c("G", str_pad(groupIndex, nchar(max(groupIndex)), pad = "0"))) |>
    select(sample, genotypeGroup)

##
## Pair status. Needs the individual label on both samples.
##
labels <- manifest |> select(sample, individual)
labelsA <- labels |> rename(sampleA = sample, individualA = individual)
labelsB <- labels |> rename(sampleB = sample, individualB = individual)
pairs <- pairs |>
    left_join(labelsA, by = "sampleA") |>
    left_join(labelsB, by = "sampleB") |>
    mutate(
        sameIndividual = individualA == individualB,
        status = case_when(
            is.na(sameIndividual) ~ NA_character_,
            is.na(identity) ~ "LOW_MARKERS",
            sameIndividual & identity < matchPct ~ "EXPECTED_MATCH_FAIL",
            !sameIndividual & identity >= matchPct ~ "UNEXPECTED_MATCH",
            !sameIndividual & identity >= relatedPct ~ "PARTIAL_MATCH",
            .default = "OK"
        )
    )

QC_LEVELS <- c("PASS", "INCOMPLETE", "WARN", "FAIL")
STATUS_LEVEL <- c(EXPECTED_MATCH_FAIL = "FAIL", UNEXPECTED_MATCH = "FAIL", PARTIAL_MATCH = "WARN")

worstOf <- function(x) {
    QC_LEVELS[max(c(1, match(x, QC_LEVELS)), na.rm = TRUE)]
}

contaminationLevelOf <- function(x) {
    case_when(
        is.na(x) ~ "INCOMPLETE",
        x >= contamFail ~ "FAIL",
        x >= contamWarn ~ "WARN",
        .default = "PASS"
    )
}

##
## Per-sample identity flags, from both sides of every flagged pair.
##
flaggedPairs <- pairs |> filter(status %in% names(STATUS_LEVEL))
pairFlags <- bind_rows(
    flaggedPairs |> transmute(sample = sampleA, other = sampleB, status, identity),
    flaggedPairs |> transmute(sample = sampleB, other = sampleA, status, identity)
) |>
    group_by(sample) |>
    summarize(
        identityFlags = str_c(glue("{status} {other} ({identity}%)"), collapse = "; "),
        identityLevel = worstOf(STATUS_LEVEL[status]),
        .groups = "drop"
    )

##
## Sample summary. The paired contamination estimate is preferred: the
## single estimate assumes heterozygous sites are 50/50, which aneuploidy
## and LOH break, so it overstates contamination in tumors and cell lines.
##
summary <- labels |>
    left_join(groups, by = "sample") |>
    left_join(samples, by = "sample") |>
    left_join(contam, by = "sample") |>
    left_join(pairFlags, by = "sample") |>
    group_by(genotypeGroup) |>
    mutate(groupMembers = map_chr(sample, \(s) str_c(setdiff(sample, s), collapse = ","))) |>
    ungroup() |>
    mutate(
        contamination = coalesce(contaminationPaired, contaminationSingle),
        contaminationMethod = case_when(
            !is.na(contaminationPaired) ~ "paired",
            !is.na(contaminationSingle) ~ "single"
        ),
        contaminationLevel = contaminationLevelOf(contamination),
        markersLevel = if_else(markersCovered < minMarkers, "WARN", "PASS"),
        identityLevel = replace_na(identityLevel, "PASS"),
        qc = pmap_chr(list(contaminationLevel, markersLevel, identityLevel),
                           \(...) worstOf(c(...))),
        contaminationFlag = if_else(
            contaminationLevel %in% c("WARN", "FAIL"),
            glue("CONTAMINATION {contamination}% ({contaminationMethod})"),
            NA_character_
        ),
        markersFlag = if_else(markersLevel == "WARN",
                              glue("LOW_MARKERS {markersCovered}"), NA_character_),
        flags = pmap_chr(list(identityFlags, contaminationFlag, markersFlag),
                         \(...) str_c(na.omit(c(...)), collapse = "; ")),
        flags = na_if(flags, "")
    ) |>
    arrange(desc(match(qc, QC_LEVELS)), match(sample, sampleOrder)) |>
    select(qc, sample, individual, genotypeGroup, groupMembers,
           contamination, contaminationMethod, partner, partnerConcordance,
           contaminationSingle, contaminationPaired,
           markersCovered, medianCoverage, pctHet, flags)

write_csv(summary, path(outDir, "conpair_summary.csv"), na = "")

pairsOut <- pairs |>
    mutate(level = STATUS_LEVEL[status]) |>
    arrange(desc(match(level, QC_LEVELS)), desc(identity)) |>
    select(sampleA, sampleB, individualA, individualB, status, identity,
           concordanceAB, nMarkersAB, concordanceBA, nMarkersBA, asymmetry,
           concordanceAll, nMarkersAll)
write_csv(pairsOut, path(outDir, "conpair_pairs.csv"), na = "")

##
## Workbook: every conpair_*.csv table as a sheet, the results first. A
## table longer than an Excel sheet allows (the concordance table, from
## about 1,024 samples) is left out and named on a note sheet. The window
## size is set so Excel does not open the workbook full screen.
##
EXCEL_MAX_ROWS <- 1048575
sheets <- list(summary = summary, pairs = pairsOut, contamination = contam,
               samples = samples, concordance = concordance, run = runTable)
tooLong <- names(keep(sheets, \(df) nrow(df) > EXCEL_MAX_ROWS))
sheets <- discard(sheets, \(df) nrow(df) > EXCEL_MAX_ROWS)
if (length(tooLong) > 0) {
    sheets$note <- tibble(note = glue("Left out, more rows than Excel allows: see conpair_{tooLong}.csv"))
}
workbook <- reduce2(
    sheets, names(sheets),
    \(wb, df, name) {
        wb |>
            openxlsx2::wb_add_worksheet(name) |>
            openxlsx2::wb_add_data(name, df, na.strings = "") |>
            openxlsx2::wb_add_font(name, dims = openxlsx2::wb_dims(x = df, select = "col_names"),
                                   bold = TRUE) |>
            openxlsx2::wb_freeze_pane(name, first_row = TRUE) |>
            openxlsx2::wb_set_col_widths(name, cols = seq_len(ncol(df)), widths = "auto")
    },
    .init = openxlsx2::wb_workbook()
) |>
    openxlsx2::wb_set_bookview(x_window = 4000, y_window = 3000,
                               window_width = 32000, window_height = 20000)
openxlsx2::wb_save(workbook, path(outDir, "conpair_tables.xlsx"), overwrite = TRUE)

##
## Text report.
##
options(width = 250)
asTable <- function(df) {
    df |>
        mutate(across(everything(), \(x) replace_na(as.character(x), ""))) |>
        as.data.frame() |>
        print(row.names = FALSE, right = FALSE) |>
        capture.output()
}

qcCounts <- tibble(qc = rev(QC_LEVELS)) |>
    left_join(count(summary, qc), by = "qc") |>
    mutate(n = replace_na(n, 0L)) |>
    glue_data("{qc} {n}") |>
    str_c(collapse = ", ")

groupLines <- summary |>
    group_by(genotypeGroup) |>
    summarize(
        samples = str_c(sample, collapse = ", "),
        individuals = str_c(sort(unique(na.omit(individual))), collapse = ", "),
        .groups = "drop"
    ) |>
    arrange(genotypeGroup)

splitIndividuals <- summary |>
    filter(!is.na(individual)) |>
    group_by(individual) |>
    filter(n_distinct(genotypeGroup) > 1) |>
    summarize(groups = str_c(sort(unique(genotypeGroup)), collapse = ", "), .groups = "drop")

mixedGroups <- summary |>
    filter(!is.na(individual)) |>
    group_by(genotypeGroup) |>
    filter(n_distinct(individual) > 1) |>
    summarize(individuals = str_c(sort(unique(individual)), collapse = ", "), .groups = "drop")

report <- c(
    glue("conpairReport.R {VERSION}   {format(Sys.time(), '%Y-%m-%d %H:%M')}"),
    glue("Output: {path_real(outDir)}"),
    glue("Samples: {length(sampleOrder)}   Markers: {run[['nMarkers']]}   ",
         "Conpair driver: {run[['version']]}"),
    glue("Genome build: {str_c(unique(manifest$build), collapse = ', ')}"),
    glue("Thresholds: match {matchPct}%, related {relatedPct}%, ",
         "contamination warn {contamWarn}% / fail {contamFail}%, ",
         "min markers {minMarkers}"),
    "",
    glue("QC: {qcCounts}"),
    if (!hasIndividual) {
        "No individual labels in the manifest: samples are grouped, not flagged."
    },
    "",
    "Flagged samples",
    if (any(summary$qc != "PASS")) {
        summary |>
            filter(qc != "PASS") |>
            select(qc, sample, individual, genotypeGroup, flags) |>
            asTable()
    } else {
        "  none"
    },
    "",
    "Genotype groups (samples whose identity is at or above the match threshold)",
    asTable(groupLines),
    "",
    "Individuals split across genotype groups",
    if (nrow(splitIndividuals) > 0) asTable(splitIndividuals) else "  none",
    "",
    "Genotype groups holding more than one individual",
    if (nrow(mixedGroups) > 0) asTable(mixedGroups) else "  none",
    "",
    "Contamination (%; paired preferred, single when no partner)",
    summary |>
        arrange(desc(contamination)) |>
        select(sample, contamination, contaminationMethod, partner,
               contaminationSingle, contaminationPaired) |>
        asTable(),
    ""
)
writeLines(report, path(outDir, "conpair_report.txt"))
cat(report, sep = "\n")

##
## Figures. One PDF, US letter landscape for every cohort size:
##   1. overview: identity of every pair, contamination of every sample;
##   2. the flagged samples as a table, left out when none is flagged;
##   3. the full concordance matrix;
##   4. one small matrix per individual, TILE_COLUMNS x TILE_ROWS a page;
##   5. contamination per sample, two columns of CONTAM_ROWS a page;
##   6. a glossary of the terms used.
## The full matrix carries sample names and values only up to
## LABEL_MAX_N samples; above that the small matrices carry them.
##
PAGE_WIDTH <- 11
PAGE_HEIGHT <- 8.5
LABEL_MAX_N <- 30
TILE_COLUMNS <- 4
TILE_ROWS <- 3
TILE_MAX_N <- 6
CONTAM_ROWS <- 50

##
## Matrix order: average-linkage clustering on 100 - identity, then
## gathered so each genotype group is contiguous.
##
nSamples <- length(sampleOrder)
distMat <- matrix(100, nSamples, nSamples, dimnames = list(sampleOrder, sampleOrder))
diag(distMat) <- 0
idPairs <- pairs |> filter(!is.na(identity))
distMat[cbind(idPairs$sampleA, idPairs$sampleB)] <- 100 - idPairs$identity
distMat[cbind(idPairs$sampleB, idPairs$sampleA)] <- 100 - idPairs$identity
hc <- hclust(as.dist(distMat), method = "average")

plotOrder <- tibble(sample = sampleOrder[hc$order], pos = seq_len(nSamples)) |>
    left_join(groups, by = "sample") |>
    group_by(genotypeGroup) |>
    mutate(groupPos = min(pos)) |>
    ungroup() |>
    arrange(groupPos, pos)
ord <- plotOrder$sample
groupSplit <- factor(plotOrder$genotypeGroup, levels = unique(plotOrder$genotypeGroup))

individualColors <- if (hasIndividual) catColors(manifest$individual) else NULL
individualOf <- deframe(labels)[ord]

##
## Concordance colors: the sequential diagnostic palette with most of its
## steps above 80%, where same-individual, partial and related pairs
## separate. Unrelated pairs (about 25-40%) take the lightest step.
##
CONCORDANCE_BREAKS <- c(40, 60, 80, 90, 95, 98, 100)
concordanceColors <- circlize::colorRamp2(CONCORDANCE_BREAKS, PAL_QC)

concordanceFill <- function(value) {
    if_else(is.na(value), NA_GREY, concordanceColors(replace_na(value, 0)))
}

## Text on a concordance cell: white on the dark steps.
concordanceInk <- function(value) {
    if_else(replace_na(value, 0) >= 90, "white", "black")
}

##
## Contamination axis: pseudo-log, so 0-1%, where most samples are, gets
## about as much room as 1-25%.
##
contamTransform <- scales::transform_pseudo_log(sigma = 0.05)
contamMax <- 1.1 * max(c(contamFail * 1.4, summary$contamination, summary$contaminationSingle),
                       na.rm = TRUE)
contamBreaks <- keep(c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100), \(b) b <= contamMax)

contamScale <- function() {
    scale_x_continuous(
        transform = contamTransform,
        breaks = contamBreaks,
        limits = c(0, contamMax),
        expand = expansion(mult = c(0.01, 0.02))
    )
}

## Shaded WARN and FAIL zones, with dashed lines at the thresholds.
contamZones <- function() {
    zones <- tibble(xmin = c(contamWarn, contamFail), xmax = c(contamFail, Inf),
                    fill = unname(PAL_PASSFAIL[c("WARN", "FAIL")]))
    list(
        geom_rect(data = zones, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                  alpha = 0.15, inherit.aes = FALSE),
        scale_fill_identity(),
        geom_vline(xintercept = c(contamWarn, contamFail), linetype = "dashed",
                   linewidth = 0.3, color = "grey40")
    )
}

formatContam <- function(x) {
    if_else(x < 10, sprintf("%.2f", x), sprintf("%.1f", x))
}

themeReport <- theme_minimal(base_size = 9) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
        plot.title.position = "plot",
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey30")
    )

pageTitleTheme <- theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    plot.caption = element_text(size = 8, color = "grey30", hjust = 0)
)

qcOf <- summary |> select(sample, qc) |> deframe()

##
## Page 1: overview.
##
## Identity of every pair, one histogram row per kind of pair. The shaded
## areas are where a pair of that kind is flagged.
##
PAIR_CLASSES <- c("same individual", "different individuals", "unlabeled", "all pairs")
## Bar color by kind of pair. Red, orange and gold are left out so that no
## bar looks like a WARN or FAIL area.
PAIR_CLASS_COLORS <- c(`same individual` = PAL_CAT[1], `different individuals` = PAL_CAT[3],
                       unlabeled = NA_GREY, `all pairs` = PAL_CAT[1])
identityPairs <- pairs |>
    mutate(
        class = case_when(
            !hasIndividual ~ "all pairs",
            is.na(sameIndividual) ~ "unlabeled",
            sameIndividual ~ "same individual",
            .default = "different individuals"
        ),
        class = fct(class, levels = PAIR_CLASSES)
    ) |>
    add_count(class, name = "nClass") |>
    mutate(panel = fct_reorder(as.character(glue("{class} ({scales::comma(nClass)} pairs)")),
                               as.integer(class)))

identityBins <- identityPairs |>
    filter(!is.na(identity)) |>
    mutate(bin = pmin(floor(identity), 99) + 0.5, class = as.character(class)) |>
    count(class, panel, bin)

identityZones <- tribble(
    ~class, ~xmin, ~xmax, ~level,
    "same individual", -Inf, matchPct, "FAIL",
    "different individuals", relatedPct, matchPct, "WARN",
    "different individuals", matchPct, Inf, "FAIL"
) |>
    inner_join(identityPairs |> distinct(class, panel) |> mutate(class = as.character(class)),
               by = "class")

LABELED_PAIRS <- 10
identityFlagged <- identityPairs |>
    filter(status %in% names(STATUS_LEVEL)) |>
    mutate(level = unname(STATUS_LEVEL[status])) |>
    arrange(desc(match(level, QC_LEVELS)), desc(identity)) |>
    mutate(label = if_else(row_number() <= LABELED_PAIRS, str_c(sampleA, " / ", sampleB), "")) |>
    left_join(identityBins |> group_by(panel) |> summarize(panelMax = max(n)), by = "panel") |>
    mutate(panelMax = replace_na(panelMax, 1))
nFlaggedPairs <- nrow(identityFlagged)

identityMin <- floor(min(c(20, identityPairs$identity), na.rm = TRUE) / 10) * 10
nLowMarkerPairs <- sum(is.na(identityPairs$identity))

flaggedPairsText <- if (!hasIndividual) {
    "No individual labels, so no pair is flagged."
} else if (nFlaggedPairs == 0) {
    "No flagged pairs."
} else if (nFlaggedPairs <= LABELED_PAIRS) {
    glue("{nFlaggedPairs} flagged pairs, each labeled at a triangle.")
} else {
    glue("{nFlaggedPairs} flagged pairs; the {LABELED_PAIRS} worst are labeled, ",
         "all are in conpair_pairs.csv.")
}
identitySubtitle <- str_c(
    "Identity: the larger of the pair's two homozygous-marker concordances. ",
    glue("Dashed lines: related {relatedPct}%, match {matchPct}%. "),
    if (hasIndividual) "A pair in a shaded area is flagged.\n" else "\n",
    flaggedPairsText,
    if (nLowMarkerPairs > 0) {
        glue(" {nLowMarkerPairs} pairs with fewer than {minMarkers} markers are not shown.")
    } else {
        ""
    }
)

identityPlot <- ggplot(identityBins, aes(bin, n)) +
    geom_rect(data = identityZones,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = level),
              alpha = 0.15, inherit.aes = FALSE) +
    geom_col(aes(fill = class), width = 0.8, show.legend = FALSE) +
    geom_vline(xintercept = c(relatedPct, matchPct), linetype = "dashed",
               linewidth = 0.3, color = "grey40") +
    geom_point(data = identityFlagged, aes(identity, 0), inherit.aes = FALSE,
               shape = 25, fill = "grey15", color = "white", size = 2.5, stroke = 0.3) +
    ggrepel::geom_text_repel(data = identityFlagged, aes(identity, 0, label = label),
                             inherit.aes = FALSE, size = 2.3, min.segment.length = 0,
                             nudge_y = 0.6 * identityFlagged$panelMax, direction = "x",
                             max.overlaps = Inf, seed = 1) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c(PAL_PASSFAIL[c("WARN", "FAIL")], PAIR_CLASS_COLORS),
                      breaks = c("WARN", "FAIL"), name = "flag") +
    scale_x_continuous(breaks = seq(0, 100, 10)) +
    coord_cartesian(xlim = c(identityMin, 100)) +
    labs(x = "Identity (%)", y = "Pairs", title = "Identity of every sample pair",
         subtitle = identitySubtitle) +
    themeReport +
    theme(legend.position = "right", strip.text = element_text(hjust = 0, face = "bold"))

##
## Contamination of every sample as a strip, one row per estimate type,
## points spread vertically so they do not overlap.
##
stripData <- summary |>
    filter(!is.na(contamination)) |>
    add_count(contaminationMethod, name = "nMethod") |>
    mutate(
        rowIndex = match(contaminationMethod, c("single", "paired")),
        rowLabel = glue("{contaminationMethod} ({nMethod})"),
        y = rowIndex + vipor::offsetX(contamTransform$transform(contamination),
                                      x = rowIndex, width = 0.32),
        label = if_else(contamination >= contamWarn, sample, "")
    )
stripRows <- stripData |> distinct(rowIndex, rowLabel)
nNoEstimate <- sum(is.na(summary$contamination))

stripPlot <- ggplot(stripData, aes(contamination, y)) +
    contamZones() +
    geom_point(shape = 21, fill = "grey30", color = "white", size = 2, stroke = 0.3) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.3, min.segment.length = 0,
                             max.overlaps = Inf, seed = 1) +
    contamScale() +
    scale_y_continuous(breaks = stripRows$rowIndex, labels = stripRows$rowLabel,
                       limits = range(stripRows$rowIndex) + c(-0.5, 0.5), expand = c(0, 0)) +
    labs(
        x = "Contamination (%)", y = NULL, title = "Contamination of every sample",
        subtitle = str_c(
            "The value used for QC: the paired estimate when the sample has a partner, ",
            "else the single estimate.\n",
            glue("Shaded: WARN from {contamWarn}%, FAIL from {contamFail}%; samples there are named."),
            if (nNoEstimate > 0) glue(" {nNoEstimate} samples with no estimate are not shown.") else ""
        )
    ) +
    themeReport +
    theme(panel.grid.major.y = element_blank())

nFlaggedSamples <- sum(summary$qc != "PASS")

groupSizes <- count(groups, genotypeGroup)
individualText <- if (hasIndividual) {
    glue("{n_distinct(na.omit(manifest$individual))} individuals")
} else {
    "no individual labels"
}
overviewSubtitle <- str_c(
    glue("{nSamples} samples, {individualText}, {nrow(groupSizes)} genotype groups ",
         "({sum(groupSizes$n > 1)} with more than one sample, {sum(groupSizes$n == 1)} with one)"),
    str_c(glue("QC: {qcCounts}"),
          if (nFlaggedSamples > 0) "; the flagged samples are listed on the next page" else ""),
    glue("Thresholds: match {matchPct}%, related {relatedPct}%, ",
         "contamination warn {contamWarn}% / fail {contamFail}%, min markers {minMarkers}. ",
         "The last page defines the terms used."),
    glue("conpairReport.R {VERSION}, {format(Sys.time(), '%Y-%m-%d %H:%M')}; ",
         "Conpair driver {run[['version']]}; build {str_c(unique(manifest$build), collapse = ', ')}"),
    path_real(outDir),
    sep = "\n"
)

overviewPage <- identityPlot / stripPlot +
    plot_layout(heights = c(2.6, 1)) +
    plot_annotation(title = "Conpair cohort report", subtitle = overviewSubtitle,
                    theme = pageTitleTheme)

##
## Page 2: the flagged samples as a table, worst first, one line per flag.
## Rows go onto pages by their estimated height; the table is built when
## the PDF device is open, so its width can be measured and its font made
## smaller when it would not fit across the page.
##
TABLE_FONT <- 12
TABLE_PAD_MM <- c(3, 2.4)
TABLE_HEIGHT_IN <- 6
TABLE_WIDTH_IN <- 10.3
TABLE_MAX_FLAG_LINES <- 12

capLines <- function(x) {
    if (length(x) <= TABLE_MAX_FLAG_LINES) {
        x
    } else {
        c(head(x, TABLE_MAX_FLAG_LINES - 1),
          str_c("... ", length(x) - TABLE_MAX_FLAG_LINES + 1, " more in conpair_summary.csv"))
    }
}

## Bin of each item when the items are packed in order into bins of the
## given capacity: rows onto pages, words onto lines. An item that does not
## fit starts the next bin, and the separator before it is dropped there;
## an item larger than a bin gets a bin of its own.
fillBins <- function(sizes, capacity, seps = 0) {
    seps <- rep_len(seps, length(sizes))
    accumulate(seq_along(sizes), \(acc, k) {
        if (acc[2] > 0 && acc[2] + seps[k] + sizes[k] > capacity) {
            c(acc[1] + 1, sizes[k])
        } else {
            c(acc[1], acc[2] + (acc[2] > 0) * seps[k] + sizes[k])
        }
    }, .init = c(1, 0))[-1] |>
        map_dbl(1)
}

tableLineIn <- TABLE_FONT * 1.2 / 72
tablePadIn <- TABLE_PAD_MM[2] / 25.4
tableRows <- summary |>
    filter(qc != "PASS") |>
    mutate(
        lines = map2(flags, contamination, \(f, x) {
            capLines(c(if (!is.na(f)) str_split_1(f, "; "),
                       if (is.na(x)) "no contamination estimate"))
        }),
        flagText = map_chr(lines, \(l) str_c(l, collapse = "\n")),
        page = fillBins(lengths(lines) * tableLineIn + tablePadIn,
                        TABLE_HEIGHT_IN - tableLineIn - tablePadIn)
    )
nTablePages <- n_distinct(tableRows$page)

flagTableSubtitle <- str_c(
    "Samples whose QC is not PASS, worst first; one line per flag, as in conpair_summary.csv.\n",
    if (hasIndividual) {
        str_c(
            glue("UNEXPECTED_MATCH: identity at or above {matchPct}% with a sample labeled as ",
                 "another individual (FAIL)."),
            "\n",
            glue("EXPECTED_MATCH_FAIL: identity below {matchPct}% with a sample labeled as ",
                 "the same individual (FAIL)."),
            "\n",
            glue("PARTIAL_MATCH: identity {relatedPct}-{matchPct}% with a sample labeled as ",
                 "another individual (WARN)."),
            "\n",
            "In these flags, the percent after the other sample is the pair's identity.\n"
        )
    },
    glue("CONTAMINATION: WARN from {contamWarn}%, FAIL from {contamFail}%; the estimate used for QC ",
         "is in parentheses."),
    "\n",
    glue("LOW_MARKERS: fewer than {minMarkers} markers covered (WARN).")
)

## Table with alternating row shading, a rule between columns and a frame.
flagTableGrob <- function(rows, fontSize) {
    df <- rows |>
        select(QC = qc, Sample = sample, Individual = individual,
               `Genotype group` = genotypeGroup, Flags = flagText) |>
        select(-any_of(if (hasIndividual) character() else "Individual")) |>
        mutate(across(everything(), \(x) replace_na(as.character(x), "")))
    nr <- nrow(df)
    nc <- ncol(df)
    tableTheme <- gridExtra::ttheme_minimal(
        base_size = fontSize,
        padding = grid::unit(TABLE_PAD_MM, "mm"),
        core = list(bg_params = list(fill = rep(rep(c("white", "grey92"), length.out = nr), nc),
                                     col = NA)),
        colhead = list(bg_params = list(fill = "grey80", col = NA),
                       fg_params = list(fontface = "bold"))
    )
    g <- gridExtra::tableGrob(df, rows = NULL, theme = tableTheme)
    ## Text left and top aligned in its cell, half the padding from the edges.
    isText <- g$layout$name %in% c("core-fg", "colhead-fg")
    g$grobs[isText] <- map(g$grobs[isText], \(t) {
        grid::editGrob(t, x = grid::unit(TABLE_PAD_MM[1] / 2, "mm"), hjust = 0,
                       y = grid::unit(1, "npc") - grid::unit(TABLE_PAD_MM[2] / 2, "mm"),
                       vjust = 1)
    })
    rule <- grid::gpar(col = "grey55", lwd = 0.7)
    g |>
        gtable::gtable_add_grob(
            map(seq_len(nc - 1), \(k) grid::segmentsGrob(1, 0, 1, 1, gp = rule)),
            t = 1, b = nr + 1, l = seq_len(nc - 1), clip = "off", name = "column-rule"
        ) |>
        gtable::gtable_add_grob(grid::segmentsGrob(0, 0, 1, 0, gp = rule),
                                t = 1, l = 1, r = nc, clip = "off", name = "header-rule") |>
        gtable::gtable_add_grob(grid::rectGrob(gp = grid::gpar(fill = NA, col = "grey40", lwd = 0.8)),
                                t = 1, b = nr + 1, l = 1, r = nc, clip = "off", name = "frame")
}

flagTablePage <- function(rows, pageNumber) {
    g <- flagTableGrob(rows, TABLE_FONT)
    widthIn <- grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE)
    if (widthIn > TABLE_WIDTH_IN) {
        g <- flagTableGrob(rows, floor(10 * TABLE_FONT * TABLE_WIDTH_IN / widthIn) / 10)
    }
    ## A null column and row take the rest of the page, so the table sits
    ## at the top left.
    g <- g |>
        gtable::gtable_add_cols(grid::unit(1, "null")) |>
        gtable::gtable_add_rows(grid::unit(1, "null"))
    title <- if (nTablePages > 1) {
        glue("Flagged samples ({nFlaggedSamples}), page {pageNumber} of {nTablePages}")
    } else {
        glue("Flagged samples ({nFlaggedSamples})")
    }
    wrap_elements(full = g) +
        plot_annotation(title = title, subtitle = flagTableSubtitle, theme = pageTitleTheme)
}

##
## Page 3: the full matrix of directed concordance. Above LABEL_MAX_N
## samples it shows the pattern only: no names except for samples whose
## QC is not PASS, no values, no group gaps. Flagged pairs are
## outlined in both directions.
##
labelMode <- nSamples <= LABEL_MAX_N

## The name of a sample that is not PASS carries why: (c) contamination at
## or above the warning threshold, (m) an identity flag, (c, m) both.
sampleCodes <- summary |>
    transmute(
        sample,
        contaminated = contaminationLevelOf(contamination) %in% c("WARN", "FAIL"),
        mismatched = sample %in% pairFlags$sample,
        code = case_when(
            contaminated & mismatched ~ "c, m",
            contaminated ~ "c",
            mismatched ~ "m",
            .default = ""
        ),
        label = if_else(code == "", sample, str_c(sample, " (", code, ")"))
    )
codeOf <- sampleCodes |> select(sample, code) |> deframe()
labelOf <- sampleCodes |> select(sample, label) |> deframe()
codeOrd <- unname(codeOf[ord])
labelOrd <- unname(labelOf[ord])

mat <- matrix(NA_real_, nSamples, nSamples, dimnames = list(sampleOrder, sampleOrder))
mat[cbind(directed$reference, directed$query)] <- directed$concHom
mat <- mat[ord, ord]

flagMat <- matrix(FALSE, nSamples, nSamples, dimnames = list(ord, ord))
flagMat[cbind(identityFlagged$sampleA, identityFlagged$sampleB)] <- TRUE
flagMat[cbind(identityFlagged$sampleB, identityFlagged$sampleA)] <- TRUE

qcOrd <- unname(qcOf[ord])
annoTable <- tibble(QC = qcOrd)
annoColors <- list(QC = PAL_PASSFAIL)
if (labelMode && !is.null(individualColors)) {
    annoTable <- annoTable |> mutate(individual = unname(individualOf))
    annoColors <- c(annoColors, list(individual = individualColors))
}

topAnnotation <- HeatmapAnnotation(
    df = as.data.frame(annoTable),
    col = annoColors,
    na_col = NA_GREY,
    annotation_name_gp = grid::gpar(fontsize = 7),
    simple_anno_size = grid::unit(2.5, "mm"),
    annotation_legend_param = list(QC = list(at = intersect(rev(QC_LEVELS), qcOrd)))
)
leftAnnotation <- rowAnnotation(
    df = as.data.frame(annoTable),
    col = annoColors,
    na_col = NA_GREY,
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(2.5, "mm")
)
markIndex <- which(qcOrd != "PASS")
rightAnnotation <- if (!labelMode && length(markIndex) > 0) {
    rowAnnotation(mark = anno_mark(at = markIndex, labels = labelOrd[markIndex],
                                   labels_gp = grid::gpar(fontsize = 7)))
} else {
    NULL
}

cellBorder <- if (nSamples <= LABEL_MAX_N) {
    grid::gpar(col = "white", lwd = 1)
} else if (nSamples <= 150) {
    grid::gpar(col = "white", lwd = 0.25)
} else {
    grid::gpar(col = NA)
}
valueSize <- if (nSamples <= 15) 7 else 5.5

matrixHeatmap <- Heatmap(
    mat,
    name = "concordance (%)",
    col = concordanceColors,
    na_col = NA_GREY,
    rect_gp = cellBorder,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = if (labelMode) groupSplit,
    column_split = if (labelMode) groupSplit,
    row_gap = grid::unit(1, "mm"),
    column_gap = grid::unit(1, "mm"),
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 7),
    column_title_gp = grid::gpar(fontsize = 7),
    show_row_names = labelMode,
    show_column_names = labelMode,
    row_labels = labelOrd,
    column_labels = labelOrd,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 7),
    top_annotation = topAnnotation,
    left_annotation = leftAnnotation,
    right_annotation = rightAnnotation,
    width = grid::unit(if (labelMode) 5.2 else 6.2, "in"),
    height = grid::unit(if (labelMode) 5.2 else 6.2, "in"),
    use_raster = nSamples > 200,
    heatmap_legend_param = list(
        at = rev(CONCORDANCE_BREAKS),
        labels = rev(c(glue("<= {CONCORDANCE_BREAKS[1]}"), CONCORDANCE_BREAKS[-1])),
        color_bar = "discrete"
    ),
    layer_fun = function(j, i, x, y, w, h, fill) {
        onDiag <- i == j
        if (any(onDiag)) {
            grid::grid.rect(x[onDiag], y[onDiag], w[onDiag], h[onDiag],
                            gp = grid::gpar(fill = "white", col = "white"))
        }
        value <- pindex(mat, i, j)
        show <- labelMode & !onDiag & !is.na(value) & value >= relatedPct
        if (any(show)) {
            grid::grid.text(sprintf("%.1f", value[show]), x[show], y[show],
                            gp = grid::gpar(fontsize = valueSize,
                                            col = concordanceInk(value[show])))
        }
        flagged <- pindex(flagMat, i, j)
        if (any(flagged)) {
            grid::grid.rect(x[flagged], y[flagged], w[flagged], h[flagged],
                            gp = grid::gpar(fill = NA, col = "black", lwd = 1.5))
        }
    }
)

matrixTitle <- str_c(
    "Concordance (%) at markers homozygous in the row (reference) sample\n",
    if (labelMode) {
        glue("Blocks are genotype groups; values shown at or above {relatedPct}%.")
    } else {
        glue("{nSamples} samples in genotype-group order. Names only for samples whose ",
             "QC is not PASS; the values for each individual are on the next pages.")
    },
    "\n",
    if (nFlaggedPairs > 0) "Outlined cells: flagged pairs." else "No flagged pairs.",
    if (any(codeOrd != "")) {
        str_c(
            "\n",
            glue("After a sample name: (c) contamination at or above {contamWarn}%, ",
                 "(m) an identity flag; the flagged-sample table gives the details."),
            if (any(qcOrd != "PASS" & codeOrd == "")) {
                "\nA flagged name with no code: too few markers or no contamination estimate."
            }
        )
    }
)

##
## Page 4: one small matrix (tile) per individual. A tile holds the
## connected samples where an edge is the same individual label or an
## identity at or above the related threshold, so a mislabeled sample
## shows up in one tile next to both its labeled individual and its
## genotype match. Samples that connect to nothing get no tile.
##
families <- igraph::graph_from_data_frame(
    pairs |>
        filter(replace_na(sameIndividual, FALSE) | replace_na(identity >= relatedPct, FALSE)) |>
        select(sampleA, sampleB),
    directed = FALSE,
    vertices = tibble(name = sampleOrder)
) |>
    igraph::components() |>
    pluck("membership") |>
    enframe(name = "sample", value = "family") |>
    left_join(summary |> select(sample, qc, individual, genotypeGroup), by = "sample") |>
    group_by(family) |>
    filter(n() > 1) |>
    summarize(
        members = list(sort(sample)),
        n = n(),
        worst = worstOf(qc),
        individuals = str_c(sort(unique(na.omit(individual))), collapse = ", "),
        groups = str_c(sort(unique(genotypeGroup)), collapse = ", "),
        firstGroup = min(genotypeGroup),
        .groups = "drop"
    ) |>
    arrange(desc(match(worst, QC_LEVELS)), firstGroup)

##
## The part of the names all samples in a tile share, cut back to its last
## separator, so "LS3601_01DR_P0_t1" and "..._t2a" become "t1" and "t2a".
##
sharedPrefix <- function(x) {
    chars <- str_split(x, "")
    width <- min(lengths(chars))
    firstDiff <- detect_index(seq_len(width), \(k) n_distinct(map_chr(chars, k)) > 1)
    common <- str_sub(x[1], 1, if (firstDiff == 0) width else firstDiff - 1)
    str_extract(common, "^.*[_.-]") |> replace_na("")
}

## sideIn is the tile's approximate side in inches; text scales with it.
makeTile <- function(members, worst, individuals, groups, sideIn) {
    n <- length(members)
    prefix <- sharedPrefix(members)
    short <- str_sub(members, str_length(prefix) + 1)
    title <- if (prefix != "") {
        str_remove(prefix, "[_.-]$")
    } else if (individuals != "") {
        individuals
    } else {
        groups
    }
    subtitle <- str_c(c(groups, worst, if (individuals != "" && individuals != title) individuals),
                      collapse = "   ")
    contamY <- -(n + 1.15)
    textMain <- min(3, 0.25 * sideIn * 25.4 / (n + 1.15))
    textSmall <- 0.85 * textMain

    cells <- expand_grid(reference = members, query = members) |>
        left_join(directed |> select(reference, query, concHom, concAll),
                  by = c("reference", "query")) |>
        mutate(
            x = match(query, members),
            y = -match(reference, members),
            onDiag = reference == query,
            fill = if_else(onDiag, "white", concordanceFill(concHom)),
            ink = concordanceInk(concHom),
            homLabel = if_else(is.na(concHom), "NA", sprintf("%.1f", concHom)),
            allLabel = if_else(is.na(concAll), "", sprintf("(%.1f)", concAll))
        )
    values <- cells |> filter(!onDiag)
    contamCells <- tibble(sample = members, x = seq_len(n)) |>
        left_join(summary |> select(sample, contamination, contaminationMethod), by = "sample") |>
        mutate(
            level = contaminationLevelOf(contamination),
            label = if_else(is.na(contamination), "none",
                            str_c(formatContam(contamination), " ",
                                  str_sub(contaminationMethod, 1, 1))),
            fill = if_else(level == "PASS", "white", unname(PAL_PASSFAIL[level])),
            ink = if_else(level %in% c("FAIL", "INCOMPLETE"), "white", "black")
        )
    rotate <- max(str_length(short)) > 5
    frame <- if (worst == "PASS") {
        element_blank()
    } else {
        element_rect(color = PAL_PASSFAIL[[worst]], fill = NA, linewidth = 1.2)
    }

    ggplot() +
        geom_tile(data = cells, aes(x, y, fill = fill), color = "white", linewidth = 0.5) +
        geom_text(data = values, aes(x, y + 0.14, label = homLabel, color = ink),
                  size = textMain) +
        geom_text(data = values, aes(x, y - 0.22, label = allLabel, color = ink),
                  size = textSmall) +
        geom_tile(data = contamCells, aes(x, contamY, fill = fill),
                  color = "grey85", linewidth = 0.3) +
        geom_text(data = contamCells, aes(x, contamY, label = label, color = ink),
                  size = textSmall) +
        scale_fill_identity() +
        scale_color_identity() +
        scale_x_continuous(breaks = seq_len(n), labels = short, position = "top",
                           expand = c(0, 0)) +
        scale_y_continuous(breaks = c(-seq_len(n), contamY), labels = c(short, "contam. %"),
                           expand = c(0, 0)) +
        coord_fixed(clip = "off") +
        labs(title = str_trunc(title, if (sideIn < 3) 30 else 80), subtitle = subtitle,
             x = NULL, y = NULL) +
        theme_minimal(base_size = 7) +
        theme(
            panel.grid = element_blank(),
            axis.text = element_text(size = 6.5),
            axis.text.x.top = element_text(angle = if (rotate) 90 else 0,
                                           hjust = if (rotate) 0 else 0.5, vjust = 0.5),
            plot.title = element_text(size = 8, face = "bold"),
            plot.subtitle = element_text(size = 6.5, color = "grey30"),
            plot.background = frame,
            plot.margin = margin(6, 6, 6, 6)
        )
}

tileArgs <- function(df) select(df, members, worst, individuals, groups)
gridFamilies <- families |> filter(n <= TILE_MAX_N)
soloFamilies <- families |> filter(n > TILE_MAX_N)
tilesPerPage <- TILE_COLUMNS * TILE_ROWS

gridPlots <- if (nrow(gridFamilies) > 0) {
    gridFamilies |>
        mutate(page = (row_number() - 1) %/% tilesPerPage + 1) %>%
        split(.$page) |>
        map(\(chunk) {
            pmap(tileArgs(chunk), makeTile, sideIn = 1.9) |>
                wrap_plots(ncol = TILE_COLUMNS, nrow = TILE_ROWS)
        })
} else {
    list()
}
soloPlots <- pmap(tileArgs(soloFamilies), makeTile, sideIn = 6.2) |>
    map(wrap_plots)
tilePlots <- unname(c(gridPlots, soloPlots))

tileSubtitle <- str_c(
    if (hasIndividual) {
        glue("Each tile holds the samples of one individual, with any sample at or above ",
             "{relatedPct}% identity to them.")
    } else {
        glue("Each tile holds samples connected by identity at or above {relatedPct}%.")
    },
    "\n",
    "Cell: concordance (%) at markers homozygous in the row (reference) sample; ",
    "below it, in parentheses, at all covered markers.\n",
    "Bottom row: contamination (%) used for QC; p = paired estimate, s = single. ",
    "Sample names omit the part they share, which is the tile title.\n",
    "Under the title: genotype group(s) and the worst QC in the tile. ",
    "A frame in a QC color marks a tile whose worst QC is not PASS."
)
tileCaption <- glue("{nSamples - sum(families$n)} samples have no tile: no other sample ",
                    "shares their individual or reaches {relatedPct}% identity with them.")
tilePages <- imap(tilePlots, \(p, k) {
    p + plot_annotation(
        title = glue("Samples of the same individual, page {k} of {length(tilePlots)}"),
        subtitle = tileSubtitle,
        caption = tileCaption,
        theme = pageTitleTheme
    )
})

##
## Page 5: contamination per sample, highest first, two columns of
## CONTAM_ROWS samples a page. Samples with no estimate come first.
## When a sample has a paired estimate as well as a single one, QC
## uses the paired one; the single one is still drawn, as a grey x, so the
## two can be compared.
##
contamRows <- summary |>
    arrange(!is.na(contamination), desc(contamination), match(sample, sampleOrder)) |>
    mutate(
        column = (row_number() - 1) %/% CONTAM_ROWS + 1,
        page = (column - 1) %/% 2 + 1,
        estimate = contaminationMethod,
        comparisonSingle = if_else(contaminationMethod == "paired", contaminationSingle, NA_real_)
    )
nContamPages <- max(contamRows$page)

ESTIMATE_SHAPES <- c(paired = 16, single = 21, comparison = 4)
ESTIMATE_COLORS <- c(paired = "grey10", single = "grey10", comparison = "grey55")
ESTIMATE_SIZES <- c(paired = 2, single = 2, comparison = 1.6)
ESTIMATE_LABELS <- c(
    paired = "Paired estimate. Used for QC.",
    single = "Single estimate, for a sample with no paired estimate. Used for QC.",
    comparison = str_c("Single estimate, for a sample that also has a paired estimate. ",
                       "For comparison only; QC uses the paired estimate.")
)
## Only the symbols that appear somewhere, in the same order on every page.
shownEstimates <- intersect(
    names(ESTIMATE_SHAPES),
    c(contamRows$estimate, if (any(!is.na(contamRows$comparisonSingle))) "comparison")
)

contamPage <- function(rows, pageNumber) {
    rows <- rows |> mutate(sample = fct_rev(fct_inorder(sample)))
    ggplot(rows, aes(y = sample)) +
        contamZones() +
        geom_segment(data = filter(rows, !is.na(comparisonSingle)),
                     aes(x = contamination, xend = comparisonSingle, yend = sample),
                     color = "grey75", linewidth = 0.3) +
        geom_point(data = filter(rows, !is.na(comparisonSingle)),
                   aes(x = comparisonSingle, shape = "comparison"),
                   color = ESTIMATE_COLORS[["comparison"]], size = 1.3) +
        geom_point(data = filter(rows, !is.na(contamination)),
                   aes(x = contamination, shape = estimate),
                   color = ESTIMATE_COLORS[["paired"]], fill = "white", size = 2, stroke = 0.6) +
        geom_text(data = filter(rows, is.na(contamination)),
                  aes(x = 0, label = "no estimate"), hjust = 0, size = 2.2, color = "grey30") +
        scale_shape_manual(values = ESTIMATE_SHAPES, breaks = shownEstimates,
                           limits = shownEstimates, labels = unname(ESTIMATE_LABELS[shownEstimates]),
                           name = NULL) +
        guides(shape = guide_legend(ncol = 1, override.aes = list(
            color = unname(ESTIMATE_COLORS[shownEstimates]),
            size = unname(ESTIMATE_SIZES[shownEstimates])
        ))) +
        contamScale() +
        facet_wrap(~column, ncol = n_distinct(rows$column), scales = "free_y") +
        labs(
            x = "Contamination (%)", y = NULL,
            title = glue("Contamination per sample, page {pageNumber} of {nContamPages}"),
            subtitle = str_c(
                glue("Highest first. Shaded: WARN from {contamWarn}%, FAIL from {contamFail}%. "),
                "A grey line joins the two estimates of a sample that has both.\n",
                "Single estimate: from the sample alone; aneuploidy and LOH inflate it. ",
                "Paired estimate: at the markers homozygous in the sample's partner, another ",
                "sample of the same individual."
            )
        ) +
        themeReport +
        theme(
            strip.text = element_blank(),
            legend.position = "top",
            legend.justification = "left",
            legend.location = "plot",
            legend.text = element_text(size = 8, margin = margin(l = 4)),
            legend.key.height = grid::unit(10, "pt"),
            legend.key.spacing.y = grid::unit(1, "pt"),
            legend.margin = margin(0, 0, 0, 0),
            axis.text.y = element_text(size = 6.5),
            panel.grid.major.y = element_line(linewidth = 0.15, color = "grey93")
        )
}

contamPages <- contamRows %>%
    split(.$page) |>
    unname() |>
    imap(contamPage)

##
## Page 6: a glossary of the terms in this report. The entries flow down
## three columns of balanced height, so no line runs across the page; an
## entry is never split, and a section heading always has its first entry
## under it. The page is built when the PDF device is open, so the text is
## wrapped by its measured width, and set smaller if it would not fit.
##
GLOSSARY_FONT <- 9
GLOSSARY_MIN_FONT <- 7
GLOSSARY_COLUMNS <- 3
GLOSSARY_WIDTH_IN <- 10.3
GLOSSARY_HEIGHT_IN <- 7.2
GLOSSARY_GUTTER_IN <- 0.35
GLOSSARY_KEY_IN <- 0.3
KEY_SIDE_IN <- 0.12

## Keys: small copies of the marks the report draws.
keySquare <- function(fill, color = NA, lwd = 1) {
    grid::rectGrob(width = grid::unit(KEY_SIDE_IN, "in"), height = grid::unit(KEY_SIDE_IN, "in"),
                   gp = grid::gpar(fill = fill, col = color, lwd = lwd))
}

## Squares side by side, one per color.
keyStrip <- function(fills, alpha = 1) {
    n <- length(fills)
    width <- 1.6 * KEY_SIDE_IN / n
    grid::rectGrob(x = grid::unit(0.5, "npc") + grid::unit((seq_len(n) - (n + 1) / 2) * width, "in"),
                   width = grid::unit(width, "in"), height = grid::unit(KEY_SIDE_IN, "in"),
                   gp = grid::gpar(fill = scales::alpha(unname(fills), alpha), col = NA))
}

## A point as ggplot2 draws it with this shape, size and stroke.
keyPoint <- function(shape, size, color, fill = NA, stroke = 0.5) {
    grid::pointsGrob(0.5, 0.5, pch = shape, size = grid::unit(1, "char"),
                     gp = grid::gpar(col = color, fill = fill,
                                     fontsize = size * .pt + stroke * .stroke / 2,
                                     lwd = stroke * .stroke / 2))
}

GLOSSARY_KEYS <- list(
    FAIL = keySquare(PAL_PASSFAIL[["FAIL"]]),
    WARN = keySquare(PAL_PASSFAIL[["WARN"]]),
    INCOMPLETE = keySquare(PAL_PASSFAIL[["INCOMPLETE"]]),
    PASS = keySquare(PAL_PASSFAIL[["PASS"]]),
    sameBar = keySquare(PAIR_CLASS_COLORS[["same individual"]]),
    otherBar = keySquare(PAIR_CLASS_COLORS[["different individuals"]]),
    unlabeledBar = keySquare(PAIR_CLASS_COLORS[["unlabeled"]]),
    zones = keyStrip(PAL_PASSFAIL[c("WARN", "FAIL")], alpha = 0.15),
    threshold = grid::segmentsGrob(0.15, 0.5, 0.85, 0.5,
                                   gp = grid::gpar(col = "grey40", lty = "dashed", lwd = 0.8)),
    flaggedPair = keyPoint(25, 2.5, "white", "grey15", stroke = 0.3),
    qcBar = keyStrip(PAL_PASSFAIL),
    outlined = keySquare(concordanceFill(70), color = "black", lwd = 1.5),
    noValue = keySquare(NA_GREY),
    paired = keyPoint(ESTIMATE_SHAPES[["paired"]], ESTIMATE_SIZES[["paired"]],
                      ESTIMATE_COLORS[["paired"]], stroke = 0.6),
    single = keyPoint(ESTIMATE_SHAPES[["single"]], ESTIMATE_SIZES[["single"]],
                      ESTIMATE_COLORS[["single"]], fill = "white", stroke = 0.6),
    comparison = keyPoint(ESTIMATE_SHAPES[["comparison"]], ESTIMATE_SIZES[["comparison"]],
                          ESTIMATE_COLORS[["comparison"]])
)

## A term of NA makes a note: the text alone, in italics.
nMarkersText <- scales::comma(as.numeric(run[["nMarkers"]]))
glossaryEntries <- tribble(
    ~section, ~term, ~text, ~key,
    "QC levels", "QC",
    "A sample's level: the worst of its identity flags, its contamination and its covered markers.", NA,
    "QC levels", "FAIL",
    glue("An UNEXPECTED_MATCH or EXPECTED_MATCH_FAIL, or contamination of {contamFail}% or more."),
    "FAIL",
    "QC levels", "WARN",
    glue("A PARTIAL_MATCH, contamination of {contamWarn}-{contamFail}%, or fewer than ",
         "{minMarkers} markers covered."),
    "WARN",
    "QC levels", "INCOMPLETE",
    "No contamination estimate could be made, and nothing is WARN or FAIL.", "INCOMPLETE",
    "QC levels", "PASS", "None of the above.", "PASS",
    "QC levels", "Flagged sample",
    "A sample whose QC is not PASS. All are listed on the Flagged samples page.", NA,

    "Flags", NA,
    "Identity flags need individual labels on both samples of a pair, and are set on both samples.",
    NA,
    "Flags", "UNEXPECTED_MATCH",
    glue("FAIL. Labeled as different individuals, identity {matchPct}% or more. Usually a sample ",
         "swap, a labeling error or one person under two IDs."),
    NA,
    "Flags", "EXPECTED_MATCH_FAIL",
    glue("FAIL. Labeled as the same individual, identity below {matchPct}%. Usually a sample swap ",
         "or a labeling error."),
    NA,
    "Flags", "PARTIAL_MATCH",
    glue("WARN. Labeled as different individuals, identity {relatedPct}-{matchPct}%. Can come from ",
         "contamination, a mixed sample, low coverage or close relatives."),
    NA,
    "Flags", "Flagged pair", "A pair with one of these three flags.", NA,
    "Flags", "CONTAMINATION",
    glue("WARN from {contamWarn}%, FAIL from {contamFail}%, by the estimate used for QC."), NA,
    "Flags", "LOW_MARKERS", glue("WARN. Fewer than {minMarkers} markers covered in the sample."), NA,
    "Flags", "(c), (m)",
    glue("After a sample name in the matrix: (c) contamination of {contamWarn}% or more, ",
         "(m) an identity flag."),
    NA,

    "Genotype and identity", "Marker",
    glue("One of the {nMarkersText} common SNPs that Conpair genotypes. Covered in a sample when ",
         "it has at least {run[['minCov']]} reads there."),
    NA,
    "Genotype and identity", "Concordance",
    "The percent of markers, covered in both samples, at which the two samples have the same genotype.",
    NA,
    "Genotype and identity", "concordanceHom",
    str_c("Concordance at the markers homozygous in the reference sample (the matrix row), ",
          "compared with the query (the column). LOH in the query does not lower it."),
    NA,
    "Genotype and identity", "concordanceAll",
    str_c("Concordance at all markers covered in both samples. The same in both directions, ",
          "but LOH lowers it. In parentheses in the tiles."),
    NA,
    "Genotype and identity", "Identity",
    str_c("The larger of a pair's two concordanceHom values, the one LOH has not lowered. ",
          "Same individual: near 100%. Unrelated: about 25-40%. Relatives: in between."),
    NA,
    "Genotype and identity", "Match threshold",
    glue("{matchPct}%. Identity at or above it means the same individual."), NA,
    "Genotype and identity", "Related threshold",
    glue("{relatedPct}%. From here up to the match threshold, identity between different ",
         "individuals is a PARTIAL_MATCH."),
    NA,
    "Genotype and identity", "Min markers",
    glue("{minMarkers}. A concordance from fewer markers is not used."), NA,
    "Genotype and identity", "Individual",
    "The person a sample is labeled as coming from, in the manifest.", NA,
    "Genotype and identity", "Genotype group",
    str_c("G01, G02, ...: samples joined by identity at or above the match threshold, directly ",
          "or through other samples. Ideally one group per individual."),
    NA,
    "Genotype and identity", "LOH",
    str_c("Loss of heterozygosity: a tumor or cell line loses one copy of a region, so its ",
          "heterozygous markers there look homozygous."),
    NA,

    "Contamination", "Contamination",
    "The estimated percent of a sample's DNA that comes from another individual.", NA,
    "Contamination", "Single estimate",
    str_c("From the sample alone, against population allele frequencies. It assumes ",
          "heterozygous markers are 50/50, so aneuploidy and LOH inflate it."),
    NA,
    "Contamination", "Paired estimate",
    str_c("Takes the markers homozygous in the sample's partner as known: reads of the other ",
          "allele there are contamination. Copy-number changes do not affect it."),
    NA,
    "Contamination", "Partner",
    glue("The other sample whose genotype best matches, at concordanceHom {matchPct}% or more. ",
         "Chosen by genotype, not by label."),
    NA,
    "Contamination", "Used for QC",
    "The paired estimate when the sample has a partner, else the single estimate.", NA,
    "Contamination", "Aneuploidy",
    "Gains or losses of chromosomes or large parts of them. Common in tumors and cell lines.", NA,
    "Contamination", "Pseudo-log axis",
    "Gives 0-1%, where most samples are, about as much room as 1-25%. The labels are true values.",
    NA,

    "Marks on the pages", "Blue bars",
    "Pairs labeled as the same individual, on page 1. Every pair when no sample is labeled.",
    "sameBar",
    "Marks on the pages", "Green bars", "Pairs labeled as different individuals.", "otherBar",
    "Marks on the pages", "Grey bars", "Pairs with an unlabeled sample.", "unlabeledBar",
    "Marks on the pages", "Shaded area",
    "A WARN (orange) or FAIL (red) range. On page 1, where a pair of that kind is flagged.",
    "zones",
    "Marks on the pages", "Dashed line", "A threshold.", "threshold",
    "Marks on the pages", "Triangle", "A flagged pair, at its identity, on page 1.", "flaggedPair",
    "Marks on the pages", "QC bar", "The strip along the matrix: each sample's QC.", "qcBar",
    "Marks on the pages", "Outlined cell", "A flagged pair, in the matrix.", "outlined",
    "Marks on the pages", "Grey cell", glue("No value: fewer than {minMarkers} markers."), "noValue",
    "Marks on the pages", "Tile",
    glue("A small matrix of one individual's samples, plus any sample at {relatedPct}% identity ",
         "or more to them. A frame in a QC color: the tile's worst QC is not PASS."),
    NA,
    "Marks on the pages", "contam. %",
    "A tile's bottom row: the contamination used for QC; p paired, s single estimate.", NA,
    "Marks on the pages", "Filled circle", "A paired estimate, used for QC.", "paired",
    "Marks on the pages", "Open circle",
    "A single estimate, used for QC: the sample has no partner.", "single",
    "Marks on the pages", "Grey x",
    "A single estimate, for comparison only: the sample's paired estimate is used.", "comparison"
)

glossaryFonts <- function(fontSize) {
    list(
        heading = grid::gpar(fontsize = fontSize + 2, fontface = "bold"),
        term = grid::gpar(fontsize = fontSize, fontface = "bold"),
        text = grid::gpar(fontsize = fontSize, col = "grey15"),
        note = grid::gpar(fontsize = fontSize, fontface = "italic", col = "grey35")
    )
}

## Width in inches of each string, in the font gp.
textWidthIn <- function(x, gp) {
    map_dbl(x, \(s) grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = gp)), "in",
                                       valueOnly = TRUE))
}

## The lines of an entry wrapped to width inches: the bold term starts the
## first line, two spaces before the text. Returns the text of each line,
## without the term.
wrapEntry <- function(term, text, width, fonts) {
    words <- str_split_1(text, " ")
    hasTerm <- !is.na(term)
    space <- textWidthIn(" ", fonts$text)
    sizes <- c(if (hasTerm) textWidthIn(term, fonts$term), textWidthIn(words, fonts$text))
    line <- fillBins(sizes, width, c(0, if (hasTerm) 2 * space, rep(space, length(words) - 1)))
    wordLine <- if (hasTerm) line[-1] else line
    map_chr(seq_len(max(line)), \(k) str_c(words[wordLine == k], collapse = " "))
}

## Every entry wrapped at fontSize, with its column and its top and bottom
## in inches below the top of the page body. An entry that opens a section
## carries the section heading.
glossaryLayout <- function(fontSize) {
    fonts <- glossaryFonts(fontSize)
    lineIn <- fontSize * 1.25 / 72
    headLineIn <- (fontSize + 2) * 1.2 / 72
    headIn <- headLineIn + 0.1
    columnIn <- (GLOSSARY_WIDTH_IN - (GLOSSARY_COLUMNS - 1) * GLOSSARY_GUTTER_IN) / GLOSSARY_COLUMNS
    blocks <- glossaryEntries |>
        mutate(
            lines = map2(term, text, \(t, x) wrapEntry(t, x, columnIn - GLOSSARY_KEY_IN, fonts)),
            opensSection = section != lag(section, default = ""),
            height = lengths(lines) * lineIn + opensSection * headIn,
            gap = if_else(opensSection, 1.2 * lineIn, 0.4 * lineIn)
        )
    ## The lowest column height that holds every entry in the columns.
    total <- sum(blocks$height) + sum(blocks$gap[-1])
    columnHeight <- detect(total / GLOSSARY_COLUMNS + seq(0, total, by = 0.01),
                           \(h) max(fillBins(blocks$height, h, blocks$gap)) <= GLOSSARY_COLUMNS)
    blocks <- blocks |>
        mutate(column = fillBins(height, columnHeight, gap)) |>
        group_by(column) |>
        mutate(bottom = cumsum(height + if_else(row_number() == 1, 0, gap)),
               top = bottom - height) |>
        ungroup()
    list(blocks = blocks, fonts = fonts, lineIn = lineIn, headLineIn = headLineIn,
         headIn = headIn, columnIn = columnIn)
}

## The grobs of one entry; top is in inches below the top of its column.
glossaryEntryGrob <- function(section, term, key, lines, opensSection, top, layout) {
    below <- \(inches) grid::unit(1, "npc") - grid::unit(inches, "in")
    fonts <- layout$fonts
    textTop <- top + opensSection * layout$headIn
    lineY <- below(textTop + (seq_along(lines) - 0.5) * layout$lineIn)
    termIn <- if (is.na(term)) 0 else textWidthIn(term, fonts$term) + 2 * textWidthIn(" ", fonts$text)
    ruleY <- below(top + layout$headLineIn + 0.03)
    grid::gTree(children = grid::gList(
        if (opensSection) {
            grid::textGrob(section, x = 0, y = below(top + layout$headLineIn / 2), hjust = 0,
                           gp = fonts$heading)
        },
        if (opensSection) {
            grid::segmentsGrob(0, ruleY, 1, ruleY, gp = grid::gpar(col = "grey60", lwd = 0.6))
        },
        if (!is.na(key)) {
            ## The key's box stops short of the term, so a wide key never touches it.
            grid::gTree(children = grid::gList(GLOSSARY_KEYS[[key]]),
                        vp = grid::viewport(x = grid::unit((GLOSSARY_KEY_IN - 0.08) / 2, "in"),
                                            y = lineY[1],
                                            width = grid::unit(GLOSSARY_KEY_IN - 0.08, "in"),
                                            height = grid::unit(layout$lineIn, "in")))
        },
        if (!is.na(term)) {
            grid::textGrob(term, x = grid::unit(GLOSSARY_KEY_IN, "in"), y = lineY[1], hjust = 0,
                           gp = fonts$term)
        },
        grid::textGrob(lines,
                       x = grid::unit(GLOSSARY_KEY_IN + c(termIn, rep(0, length(lines) - 1)), "in"),
                       y = lineY, hjust = 0, gp = if (is.na(term)) fonts$note else fonts$text)
    ))
}

glossaryPage <- function() {
    fontSize <- detect(seq(GLOSSARY_FONT, GLOSSARY_MIN_FONT, by = -0.5),
                       \(s) max(glossaryLayout(s)$blocks$bottom) <= GLOSSARY_HEIGHT_IN,
                       .default = GLOSSARY_MIN_FONT)
    layout <- glossaryLayout(fontSize)
    columns <- layout$blocks |>
        mutate(grob = pmap(list(section, term, key, lines, opensSection, top),
                           \(...) glossaryEntryGrob(..., layout = layout))) %>%
        split(.$column) |>
        map(\(b) grid::gTree(children = rlang::exec(grid::gList, !!!b$grob)))
    nColumns <- length(columns)
    widths <- rep(c(layout$columnIn, GLOSSARY_GUTTER_IN), length.out = 2 * nColumns - 1)
    ## A null column takes the rest of the width, so the text sits at the left.
    body <- gtable::gtable(widths = grid::unit(c(widths, 1), c(rep("in", length(widths)), "null")),
                           heights = grid::unit(1, "null")) |>
        gtable::gtable_add_grob(columns, t = 1, l = seq(1, by = 2, length.out = nColumns),
                                clip = "off", name = "glossary-column")
    wrap_elements(full = body) +
        plot_annotation(title = "Glossary",
                        subtitle = "The terms used in this report, with this run's thresholds.",
                        theme = pageTitleTheme)
}

pdf(path(outDir, "conpair_report.pdf"), width = PAGE_WIDTH, height = PAGE_HEIGHT)
print(overviewPage)
tableRows %>%
    split(.$page) |>
    unname() |>
    iwalk(\(rows, k) print(flagTablePage(rows, k)))
draw(matrixHeatmap, column_title = matrixTitle,
     column_title_gp = grid::gpar(fontsize = 10), merge_legends = TRUE)
walk(tilePages, print)
walk(contamPages, print)
print(glossaryPage())
invisible(dev.off())

quit(status = 0)
