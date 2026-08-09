#!/usr/bin/env bash
set -euo pipefail

SIM_DIR="${SILK_SIM_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIM_DIR="$(cd "$SIM_DIR" && pwd)"
export SILK_SIM_DIR="$SIM_DIR"
export SILK_OUT_DIR="${SILK_OUT_DIR:-$SIM_DIR/outputs}"

SILK_PACKAGE_REPOSITORY="${SILK_PACKAGE_REPOSITORY:-https://github.com/mehrdadmhmdi/SILK.git}"
SILK_PACKAGE_REF="${SILK_PACKAGE_REF:-master}"
SILK_INSTALL_ID="${SLURM_JOB_ID:-manual-$$}"
PACKAGE_SOURCE_PARENT="$SILK_OUT_DIR/package-source/$SILK_INSTALL_ID"
PACKAGE_SOURCE="$PACKAGE_SOURCE_PARENT/SILK"
export SILK_R_LIB="${SILK_R_LIB:-$SILK_OUT_DIR/R-library/$SILK_INSTALL_ID}"

mkdir -p "$PACKAGE_SOURCE_PARENT" "$SILK_R_LIB" "$SILK_OUT_DIR/summary"

# Preserve the library paths R would use before introducing the run-specific
# package library. Cluster dependencies can therefore remain in the Cray R or
# user libraries, while SILK itself is loaded from this run's private library.
existing_r_libraries="$(Rscript -e 'cat(.libPaths(), collapse=.Platform$path.sep)')"

echo "Cloning SILK once from $SILK_PACKAGE_REPOSITORY (ref $SILK_PACKAGE_REF)"
git clone --depth 1 --single-branch --branch "$SILK_PACKAGE_REF" \
  "$SILK_PACKAGE_REPOSITORY" "$PACKAGE_SOURCE"
export SILK_PACKAGE_COMMIT="$(git -C "$PACKAGE_SOURCE" rev-parse HEAD)"

echo "Installing GitHub commit $SILK_PACKAGE_COMMIT into $SILK_R_LIB"
R CMD INSTALL --preclean --no-multiarch --library="$SILK_R_LIB" "$PACKAGE_SOURCE"
export R_LIBS_USER="$SILK_R_LIB${existing_r_libraries:+:${existing_r_libraries}}"

cat > "$SILK_OUT_DIR/summary/silk_package_provenance.txt" <<EOF
repository=$SILK_PACKAGE_REPOSITORY
ref=$SILK_PACKAGE_REF
commit=$SILK_PACKAGE_COMMIT
source=$PACKAGE_SOURCE
library=$SILK_R_LIB
EOF

# Fail before array submission if GitHub does not contain the characteristic-
# kernel package API required by these simulation scripts.
Rscript -e '
  library(SILK)
  defaults <- silk_default_options()
  stopifnot(
    is.function(fit_silk_registration),
    is.function(predict_silk_registration),
    identical(defaults$BIOMARKER_BANDWIDTH, "median"),
    !any(grepl("Linear", defaults$METHODS, fixed = TRUE))
  )
  cat(
    "Loaded SILK", as.character(packageVersion("SILK")),
    "from", find.package("SILK"),
    "with exact Gaussian biomarker-kernel registration\n"
  )
'
