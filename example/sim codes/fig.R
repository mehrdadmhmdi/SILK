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
