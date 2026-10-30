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
manifest.tsv and the conpair_*.tsv tables written by conpairCohort.py.

Options (all in percent except --min-markers):
  --match <x>        Identity at or above this means the same individual.
                     Default: the value conpairCohort.py used, from
                     conpair_run.tsv, so partners and groups agree.
  --related <x>      Identity at or above this between samples labeled as
                     different individuals is a PARTIAL_MATCH warning.
                     Default: 80
  --min-markers <n>  A concordance based on fewer markers is not used.
                     Default: the value conpairCohort.py used.
  --contam-warn <x>  Contamination warning threshold. Default: 2
  --contam-fail <x>  Contamination failure threshold. Default: 5
  -h, --help         Show this message and exit.

Writes into <OutDir>:
  conpair_summary.tsv   one row per sample, worst verdict first
  conpair_pairs.tsv     one row per sample pair, flagged pairs first
  conpair_report.txt    text report, also printed to stdout
  conpair_report.pdf    concordance heatmaps and contamination plot
") |> cat("\n")
    quit(status = 0)
}

if (length(args) == 0 || any(args %in% c("-h", "--help")) || str_starts(args[1], "-")) {
    usage()
}
outDir <- args[1]

readTable <- function(name) {
    read_tsv(path(outDir, name), na = "NA", show_col_types = FALSE)
}

run <- read_tsv(path(outDir, "conpair_run.tsv"),
                col_types = cols(.default = col_character()),
                show_col_types = FALSE) |>
    deframe()

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

samples <- readTable("conpair_samples.tsv")
contam <- readTable("conpair_contamination.tsv") |>
    mutate(partner = as.character(partner))

##
## Directed concordance. A value resting on fewer than minMarkers markers
## is dropped rather than trusted.
##
directed <- readTable("conpair_concordance.tsv") |>
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

VERDICTS <- c("PASS", "INCOMPLETE", "WARN", "FAIL")
STATUS_LEVEL <- c(EXPECTED_MATCH_FAIL = "FAIL", UNEXPECTED_MATCH = "FAIL", PARTIAL_MATCH = "WARN")

worstOf <- function(x) {
    VERDICTS[max(c(1, match(x, VERDICTS)), na.rm = TRUE)]
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
        contaminationLevel = case_when(
            is.na(contamination) ~ "INCOMPLETE",
            contamination >= contamFail ~ "FAIL",
            contamination >= contamWarn ~ "WARN",
            .default = "PASS"
        ),
        markersLevel = if_else(markersCovered < minMarkers, "WARN", "PASS"),
        identityLevel = replace_na(identityLevel, "PASS"),
        verdict = pmap_chr(list(contaminationLevel, markersLevel, identityLevel),
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
    arrange(desc(match(verdict, VERDICTS)), match(sample, sampleOrder)) |>
    select(verdict, sample, individual, genotypeGroup, groupMembers,
           contamination, contaminationMethod, partner, partnerConcordance,
           contaminationSingle, contaminationPaired,
           markersCovered, medianCoverage, pctHet, flags)

write_tsv(summary, path(outDir, "conpair_summary.tsv"), na = "")

pairsOut <- pairs |>
    mutate(level = STATUS_LEVEL[status]) |>
    arrange(desc(match(level, VERDICTS)), desc(identity)) |>
    select(sampleA, sampleB, individualA, individualB, status, identity,
           concordanceAB, nMarkersAB, concordanceBA, nMarkersBA, asymmetry,
           concordanceAll, nMarkersAll)
write_tsv(pairsOut, path(outDir, "conpair_pairs.tsv"), na = "")

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

verdictCounts <- tibble(verdict = rev(VERDICTS)) |>
    left_join(count(summary, verdict), by = "verdict") |>
    mutate(n = replace_na(n, 0L)) |>
    glue_data("{verdict} {n}") |>
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
    glue("Thresholds: match {matchPct}%, related {relatedPct}%, ",
         "contamination warn {contamWarn}% / fail {contamFail}%, ",
         "min markers {minMarkers}"),
    "",
    glue("Verdicts: {verdictCounts}"),
    if (!hasIndividual) {
        "No individual labels in the manifest: samples are grouped, not flagged."
    },
    "",
    "Flagged samples",
    if (any(summary$verdict != "PASS")) {
        summary |>
            filter(verdict != "PASS") |>
            select(verdict, sample, individual, genotypeGroup, flags) |>
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
## Figures. Samples are ordered by average-linkage clustering on
## 100 - identity, then gathered so each genotype group is contiguous.
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
nGroups <- nlevels(groupSplit)

individualColors <- if (hasIndividual) catColors(manifest$individual) else NULL
individualOf <- deframe(labels)[ord]

##
## Concordance colors: the sequential diagnostic palette with most of its
## steps above 80%, where same-individual, partial and related pairs
## separate. Unrelated pairs (about 25-40%) take the lightest step.
##
CONCORDANCE_BREAKS <- c(40, 60, 80, 90, 95, 98, 100)
concordanceColors <- circlize::colorRamp2(CONCORDANCE_BREAKS, PAL_QC)
showValues <- nSamples <= 40
fontSize <- if (nSamples <= 20) 8 else if (nSamples <= 60) 6 else 4

concordanceHeatmap <- function(valueCol, legendTitle) {
    mat <- matrix(NA_real_, nSamples, nSamples, dimnames = list(sampleOrder, sampleOrder))
    mat[cbind(directed$reference, directed$query)] <- directed[[valueCol]]
    mat <- mat[ord, ord]

    annotations <- if (!is.null(individualColors)) {
        list(
            top = HeatmapAnnotation(
                individual = individualOf,
                col = list(individual = individualColors),
                na_col = NA_GREY,
                annotation_name_gp = grid::gpar(fontsize = fontSize),
                simple_anno_size = grid::unit(3, "mm")
            ),
            left = rowAnnotation(
                individual = individualOf,
                col = list(individual = individualColors),
                na_col = NA_GREY,
                show_legend = FALSE,
                show_annotation_name = FALSE,
                simple_anno_size = grid::unit(3, "mm")
            )
        )
    } else {
        list(top = NULL, left = NULL)
    }

    Heatmap(
        mat,
        name = legendTitle,
        col = concordanceColors,
        na_col = NA_GREY,
        rect_gp = grid::gpar(col = "white", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = groupSplit,
        column_split = groupSplit,
        row_gap = grid::unit(1, "mm"),
        column_gap = grid::unit(1, "mm"),
        row_title_rot = 0,
        row_title_gp = grid::gpar(fontsize = fontSize),
        column_title_gp = grid::gpar(fontsize = fontSize),
        row_title = if (nGroups <= 30) "%s" else NULL,
        column_title = if (nGroups <= 30) "%s" else NULL,
        row_names_gp = grid::gpar(fontsize = fontSize),
        column_names_gp = grid::gpar(fontsize = fontSize),
        top_annotation = annotations$top,
        left_annotation = annotations$left,
        heatmap_legend_param = list(
            at = rev(CONCORDANCE_BREAKS),
            labels = rev(c(glue("<= {CONCORDANCE_BREAKS[1]}"), CONCORDANCE_BREAKS[-1])),
            color_bar = "discrete"
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            value <- mat[i, j]
            if (rownames(mat)[i] == colnames(mat)[j]) {
                grid::grid.rect(x, y, width, height,
                                gp = grid::gpar(fill = "white", col = "white"))
            } else if (showValues && !is.na(value) && value >= relatedPct) {
                grid::grid.text(
                    sprintf("%.1f", value), x, y,
                    gp = grid::gpar(fontsize = fontSize - 1,
                                    col = if (value >= 90) "white" else "black")
                )
            }
        }
    )
}

contamLong <- summary |>
    select(sample, contamination, paired = contaminationPaired, single = contaminationSingle) |>
    pivot_longer(c(paired, single), names_to = "estimate", values_to = "pct") |>
    filter(!is.na(pct)) |>
    mutate(sample = fct_reorder(sample, contamination, .na_rm = TRUE))
estimateColors <- catColors(contamLong$estimate)
xMax <- max(c(contamFail * 1.4, contamLong$pct * 1.05), na.rm = TRUE)

contamPlot <- ggplot(contamLong, aes(pct, sample, color = estimate, shape = estimate)) +
    geom_vline(xintercept = c(contamWarn, contamFail), linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    annotate("text", x = c(contamWarn, contamFail), y = Inf,
             label = c(glue("warn {contamWarn}%"), glue("fail {contamFail}%")),
             vjust = -0.4, hjust = 0.5, size = 3, color = "grey30") +
    geom_point(size = 2.5) +
    scale_color_manual(values = estimateColors) +
    scale_shape_manual(values = c(paired = 16, single = 1)) +
    scale_x_continuous(limits = c(0, xMax), expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(clip = "off") +
    labs(
        x = "Contamination (%)", y = NULL, color = NULL, shape = NULL,
        title = "Contamination per sample",
        subtitle = str_c(
            "paired: Conpair tumor-sample estimate, using the partner's homozygous markers\n",
            "single: Conpair normal-sample estimate from the sample alone; ",
            "overstated by aneuploidy and LOH"
        )
    ) +
    theme_minimal(base_size = 10) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
        panel.grid.major.x = element_line(linewidth = 0.2, color = "grey90"),
        legend.position = "top",
        plot.title.position = "plot",
        plot.margin = margin(18, 10, 5, 5)
    )

side <- max(7, 3 + 0.28 * nSamples)
pdf(path(outDir, "conpair_report.pdf"), width = side + 2, height = side)
draw(
    concordanceHeatmap("concHom", "concordance (%)"),
    column_title = str_c(
        "Concordance at markers homozygous in the row (reference) sample\n",
        "Blocks are genotype groups; the color bars give the expected individual"
    ),
    column_title_gp = grid::gpar(fontsize = 10)
)
print(contamPlot)
draw(
    concordanceHeatmap("concAll", "concordance (%)"),
    column_title = str_c(
        "Concordance at all covered markers (symmetric)\n",
        "LOH lowers this between samples of the same individual"
    ),
    column_title_gp = grid::gpar(fontsize = 10)
)
invisible(dev.off())

quit(status = 0)
