#!/usr/bin/env bash
set -euo pipefail

SIM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "$SIM_DIR/../.." && pwd)"
export SILK_OUT_DIR="${SILK_OUT_DIR:-$SIM_DIR/outputs}"
export SILK_R_LIB="${SILK_R_LIB:-$SILK_OUT_DIR/R-library}"
mkdir -p "$SILK_R_LIB"

# Preserve the library paths R would use before introducing the run-specific
# package library. This matters on clusters where dependencies live in the
# user's default R library even though R_LIBS_USER is not explicitly exported.
existing_r_libraries="$(Rscript -e 'cat(.libPaths(), collapse=.Platform$path.sep)')"

echo "Installing SILK from $PACKAGE_ROOT into $SILK_R_LIB"
R CMD INSTALL --preclean --no-multiarch --library="$SILK_R_LIB" "$PACKAGE_ROOT"
export R_LIBS_USER="$SILK_R_LIB${existing_r_libraries:+:${existing_r_libraries}}"

Rscript -e 'library(SILK); defaults <- silk_default_options(); stopifnot(is.function(fit_silk_registration), is.function(predict_silk_registration), identical(defaults$BIOMARKER_BANDWIDTH, "median")); cat("Installed SILK", as.character(packageVersion("SILK")), "with exact Gaussian biomarker-kernel registration\n")'
