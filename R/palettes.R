##
## Color palettes shared by every BAMQC figure. Defined once here as hex
## vectors so figures do not change when a palette package is upgraded.
## Do not hard-code colors in individual figures.
##

## Sequential, diagnostics: grey -> gold -> red -> deep purple
PAL_QC <- c("#F2F2F2", "#FFE0A3", "#FDA43A", "#F03B20",
            "#BD0026", "#7A0177", "#49006A")

## Categorical, up to 23 levels: pals::cols25() without its two greys
PAL_CAT <- c("#1F78C8", "#FF0000", "#33A02C", "#6A33C2", "#FF7F00",
             "#FFD700", "#A6CEE3", "#FB6496", "#B2DF8A", "#CAB2D6",
             "#FDBF6F", "#EEE685", "#C8308C", "#FF83FA", "#C814FA",
             "#0000FF", "#36648B", "#00E2E5", "#00FF00", "#778B00",
             "#BEBE00", "#8B3B00", "#A52A3C")

## Categorical, 24-32 levels: pals::glasbey(32). Replaces PAL_CAT
## entirely; do not append it to PAL_CAT.
PAL_CAT32 <- c("#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6",
               "#005300", "#FFD300", "#009FFF", "#9A4D42", "#00FFBE",
               "#783FC1", "#1F9698", "#FFACFD", "#B1CC71", "#F1085C",
               "#FE8F42", "#DD00FF", "#201A01", "#720055", "#766C95",
               "#02AD24", "#C8FF00", "#886C00", "#FFB79F", "#858567",
               "#A10300", "#14F9FF", "#00479E", "#DC5E93", "#93D4FF",
               "#004CFF", "#F2F318")

## "None" / unassigned
NA_GREY <- "grey70"

## QC levels. WARN and FAIL are status warning and critical colors; PASS is
## a quiet grey so problems stand out; INCOMPLETE (no estimate) is a dark
## grey. Always shown with the level's word, never by color alone.
PAL_PASSFAIL <- c(PASS = "#A6A6A6", INCOMPLETE = "#5E5E5E",
                  WARN = "#FAB219", FAIL = "#D03B3B")

##
## Named categorical colors for a vector of labels. Levels are ordered by
## descending count, ties by first appearance, so the largest groups get
## the first colors. NA is not a level; callers color it NA_GREY.
## Returns NULL past 32 levels, where no palette stays distinguishable.
##
catColors <- function(x) {
    lv <- tibble::tibble(level = x) |>
        dplyr::filter(!is.na(level)) |>
        dplyr::mutate(first = dplyr::row_number()) |>
        dplyr::group_by(level) |>
        dplyr::summarize(n = dplyr::n(), first = min(first)) |>
        dplyr::arrange(dplyr::desc(n), first) |>
        dplyr::pull(level)
    if (length(lv) <= length(PAL_CAT)) {
        rlang::set_names(PAL_CAT[seq_along(lv)], lv)
    } else if (length(lv) <= length(PAL_CAT32)) {
        rlang::set_names(PAL_CAT32[seq_along(lv)], lv)
    } else {
        NULL
    }
}
