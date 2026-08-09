#!/usr/bin/env bash
set -euo pipefail

SIM_DIR="${SILK_SIM_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIM_DIR="$(cd "$SIM_DIR" && pwd)"
export SILK_SIM_DIR="$SIM_DIR"
export SILK_OUT_DIR="${SILK_OUT_DIR:-$SIM_DIR/outputs}"

# This environment is persistent across controller and array jobs. Dependencies
# and SILK are installed only when absent or when the requested GitHub commit
# changes; array workers never install software.
export SILK_ENV_ROOT="${SILK_ENV_ROOT:-$SIM_DIR/.cluster-env-r44}"
export SILK_R_LIB="${SILK_R_LIB:-$SILK_ENV_ROOT/R-library}"
export SILK_PY_ENV="${SILK_PY_ENV:-$SILK_ENV_ROOT/pycox-cpu}"
PACKAGE_SOURCE="$SILK_ENV_ROOT/package-source/SILK"
PROVENANCE_DIR="$SILK_ENV_ROOT/provenance"
SILK_PACKAGE_REPOSITORY="${SILK_PACKAGE_REPOSITORY:-https://github.com/mehrdadmhmdi/SILK.git}"
SILK_PACKAGE_REF="${SILK_PACKAGE_REF:-master}"

mkdir -p "$SILK_R_LIB" "$PROVENANCE_DIR" "$SILK_OUT_DIR/summary"
export R_LIBS="$SILK_R_LIB${R_LIBS:+:$R_LIBS}"

if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 is required to build the persistent CPU DeepSurv environment." >&2
  exit 1
fi

if [[ ! -x "$SILK_PY_ENV/bin/python" ]]; then
  python3 -m venv "$SILK_PY_ENV"
  "$SILK_PY_ENV/bin/python" -m pip install --upgrade pip setuptools wheel
fi

if ! "$SILK_PY_ENV/bin/python" -c 'import torch' >/dev/null 2>&1; then
  "$SILK_PY_ENV/bin/python" -m pip install \
    --index-url https://download.pytorch.org/whl/cpu torch
fi
if ! "$SILK_PY_ENV/bin/python" -c 'import pycox, torchtuples, numpy, pandas, sklearn' >/dev/null 2>&1; then
  "$SILK_PY_ENV/bin/python" -m pip install \
    pycox torchtuples numpy pandas scipy scikit-learn
fi

export RETICULATE_PYTHON="$SILK_PY_ENV/bin/python"
export PYTHONNOUSERSITE=1
export CUDA_VISIBLE_DEVICES=""

# Install missing R dependencies into the persistent private library. Default
# site and user libraries remain visible, so already installed packages such as
# rlang are reused instead of hidden.
Rscript --vanilla -e '
  private <- Sys.getenv("SILK_R_LIB")
  .libPaths(c(private, .libPaths()))
  required <- c(
    "survival", "rlang", "nlme", "ggplot2", "ranger",
    "reticulate", "survivalmodels"
  )
  if (identical(tolower(Sys.getenv("SILK_ENABLE_JMBAYES2", "true")), "true")) {
    required <- c(required, "JMbayes2")
  }
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    install.packages(missing, lib = private, dependencies = NA,
                     Ncpus = max(1L, min(4L, parallel::detectCores())))
  }
  still_missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing)) {
    stop("Missing required R packages after setup: ", paste(still_missing, collapse = ", "))
  }
'

if [[ ! -d "$PACKAGE_SOURCE/.git" ]]; then
  mkdir -p "$(dirname "$PACKAGE_SOURCE")"
  git clone "$SILK_PACKAGE_REPOSITORY" "$PACKAGE_SOURCE"
fi
git -C "$PACKAGE_SOURCE" fetch --depth 1 origin "$SILK_PACKAGE_REF"
git -C "$PACKAGE_SOURCE" checkout --detach FETCH_HEAD
export SILK_PACKAGE_COMMIT="$(git -C "$PACKAGE_SOURCE" rev-parse HEAD)"

installed_commit=""
if [[ -f "$PROVENANCE_DIR/silk-installed-commit.txt" ]]; then
  installed_commit="$(<"$PROVENANCE_DIR/silk-installed-commit.txt")"
fi
if [[ "$installed_commit" != "$SILK_PACKAGE_COMMIT" ]] || \
   [[ ! -d "$SILK_R_LIB/SILK" ]]; then
  echo "Installing SILK GitHub commit $SILK_PACKAGE_COMMIT into $SILK_R_LIB"
  R CMD INSTALL --preclean --no-multiarch --library="$SILK_R_LIB" "$PACKAGE_SOURCE"
  printf '%s\n' "$SILK_PACKAGE_COMMIT" > "$PROVENANCE_DIR/silk-installed-commit.txt"
else
  echo "Reusing installed SILK GitHub commit $SILK_PACKAGE_COMMIT"
fi

cat > "$SILK_OUT_DIR/summary/silk_package_provenance.txt" <<EOF
repository=$SILK_PACKAGE_REPOSITORY
ref=$SILK_PACKAGE_REF
commit=$SILK_PACKAGE_COMMIT
source=$PACKAGE_SOURCE
library=$SILK_R_LIB
python=$RETICULATE_PYTHON
EOF

"$SILK_PY_ENV/bin/python" -m pip freeze > "$PROVENANCE_DIR/python-freeze.txt"

# Fail before array submission if either the package API or Python backend is
# incomplete. This catches configuration failures before thousands of tasks run.
Rscript --vanilla -e '
  library(SILK)
  stopifnot(
    is.function(fit_silk_registration),
    is.function(predict_silk_registration),
    identical(silk_default_options()$BIOMARKER_KERNEL, "gaussian"),
    !any(grepl("Linear", silk_default_options()$METHODS, fixed = TRUE))
  )
  stopifnot(
    reticulate::py_module_available("torch"),
    reticulate::py_module_available("pycox"),
    reticulate::py_module_available("torchtuples")
  )
  torch <- reticulate::import("torch")
  stopifnot(!isTRUE(torch$cuda$is_available()))
  survivalmodels::set_seed(20260809L)
  smoke <- data.frame(
    time = seq(1, 40),
    status = rep(c(1L, 0L), 20),
    x1 = stats::rnorm(40),
    x2 = stats::rnorm(40)
  )
  deep_fit <- survivalmodels::deepsurv(
    survival::Surv(time, status) ~ x1 + x2,
    data = smoke,
    frac = 0.2,
    num_nodes = c(4L),
    epochs = 1L,
    batch_size = 16L,
    device = "cpu",
    num_workers = 0L,
    verbose = FALSE
  )
  deep_prediction <- as.matrix(stats::predict(
    deep_fit, newdata = smoke[1:2, c("x1", "x2")], type = "survival"
  ))
  stopifnot(any(dim(deep_prediction) == 2L), all(is.finite(deep_prediction)))
  cat("Loaded SILK", as.character(packageVersion("SILK")), "from", find.package("SILK"), "\n")
  cat("DeepSurv Python:", reticulate::py_config()$python, "(CPU only)\n")
'
