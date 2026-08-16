#!/usr/bin/env Rscript
# =============================================================================
# fig.R
# Create prediction-performance plots from collected simulation summaries.
# =============================================================================

source(file.path("R", "run_simulation_comparisons.R"))
source_silk_prediction_modules()
# Use the plotting code from this checkout. The installed package may predate
# presentation-only label and palette revisions, while all statistical
# summaries remain the collected confirmatory outputs.
source(file.path("..", "..", "R", "plotting.R"))

# Comparator colors belong to the simulation/application layer.  The package
# palette intentionally contains SILK methods only.
method_palette <- function() {
  method_colours <- c(
    "SILK-Gaussian" = UIUC_ORANGE,
    "SILK-Laplace" = UIUC_ORANGE,
    "SILK-Matern32" = UIUC_ORANGE,
    "Recorded-Cox" = "#707372",
    "Recorded-Beran" = "#8E9090",
    "MMLM-Correct" = "#007E8E",
    "MMLM-Misspecified" = "#4DB6C2",
    "JM-Correct" = "#5C0E41",
    "JM-Misspecified" = "#A64D79",
    "DeepSurv" = "#006230",
    "RSF" = "#1D58A7",
    "TimeError-Integrated-Landmark" = "#7D3E13",
    "Oracle-Cox" = "#000000",
    "Oracle-Beran" = UIUC_BLUE
  )
  stats::setNames(unname(method_colours[METHOD_ORDER]), pretty_method(METHOD_ORDER))
}

method_linetype_palette <- function() {
  styles <- stats::setNames(rep("solid", length(METHOD_ORDER)), METHOD_ORDER)
  styles["SILK-Laplace"] <- "dashed"
  styles["SILK-Matern32"] <- "dotted"
  styles["Recorded-Beran"] <- "dotdash"
  styles["MMLM-Misspecified"] <- "dashed"
  styles["JM-Misspecified"] <- "dashed"
  styles["Oracle-Cox"] <- "longdash"
  styles["Oracle-Beran"] <- "twodash"
  stats::setNames(unname(styles[METHOD_ORDER]), pretty_method(METHOD_ORDER))
}

method_labels <- pretty_method(METHOD_ORDER)
if (anyDuplicated(method_labels)) {
  stop("Figure method labels are not unique.", call. = FALSE)
}
method_colours <- method_palette()
if (anyNA(method_colours)) {
  stop(
    "Figure colours are missing for: ",
    paste(names(method_colours)[is.na(method_colours)], collapse = ", "),
    call. = FALSE
  )
}

make_prediction_plots(OUT_DIR)
cat("Prediction figures written under ", file.path(OUT_DIR, "figures"), "\n", sep = "")
