#!/usr/bin/env Rscript
# =============================================================================
# fig.R
# Create prediction-performance plots from collected simulation summaries.
# =============================================================================

source(file.path("R", "run_simulation_comparisons.R"))
source_silk_prediction_modules()

make_prediction_plots(OUT_DIR)
cat("Prediction figures written under ", file.path(OUT_DIR, "figures"), "\n", sep = "")
