# =============================================================================
# dgm_biomarker_information.R
#
# How much latent-stage information do the biomarker features carry in each
# scenario of the current data-generating mechanism?
#
# Why this matters. SILK's estimand conditions on the CALIBRATED age and the
# shift-invariant history, not on the latent age; this script is not a check on
# SILK. It characterizes the DGM. The current MMLM, JM, RSF, DeepSurv, and
# time-error comparators use B1, whereas SILK uses B1--B4 jointly. It would
# therefore be incorrect to use an all-marker random forest as evidence that
# every comparator can recover stage. This script reports B1-only and all-marker
# results separately, for both landmark-only and full-history features.
#
# Reported per scenario:
#   rmse_b1_landmark       RMSE of A* from landmark B1
#   rmse_b1_full_history   RMSE of A* from B1 at all visits
#   rmse_all_landmark      RMSE of A* from landmark B1--B4
#   rmse_all_full_history  RMSE of A* from B1--B4 at all visits
#   r2_*                   corresponding out-of-sample R^2 values
#   sigma_eps           the scenario's origin-error scale, for reference
#
# Estimation uses an out-of-sample random-forest regression so the answer does
# not depend on knowing the parametric form of the biomarker mean.
#
# Interpretation limits. This is an ACHIEVABILITY diagnostic under one specific
# learner, not an information-theoretic bound: a different learner could do
# better, so the reported RMSE is an upper bound on the attainable error and the
# reported R^2 is a lower bound on the recoverable signal. Results are averaged
# over several independent train/test replications and reported with a Monte
# Carlo standard error, because a single split is noisy.
#
# Usage:
#   Rscript dgm_biomarker_information.R [n] [n_rep] [seed] [output_csv]
# =============================================================================

suppressPackageStartupMessages({
  library(SILK)
})
source(file.path("R", "cfg.R"))

args <- commandArgs(trailingOnly = TRUE)
N <- if (length(args) >= 1L) as.integer(args[[1]]) else 4000L
N_REP <- if (length(args) >= 2L) as.integer(args[[2]]) else 5L
SEED <- if (length(args) >= 3L) as.integer(args[[3]]) else 20260810L
OUT <- if (length(args) >= 4L) args[[4]] else "dgm_biomarker_information.csv"

has_ranger <- requireNamespace("ranger", quietly = TRUE)
if (!has_ranger) {
  message("ranger not available; falling back to an ordinary linear model ",
          "(main effects only). The fallback is markedly weaker for the ",
          "nonlinear biomarker means in this DGM, so treat its RMSE as a loose ",
          "upper bound and install ranger before reporting these numbers.")
}

fit_predict_rmse <- function(x_train, y_train, x_test, y_test) {
  df_tr <- data.frame(y = y_train, x_train, check.names = TRUE)
  df_te <- data.frame(x_test, check.names = TRUE)
  if (has_ranger) {
    m <- ranger::ranger(y ~ ., data = df_tr, num.trees = 500L,
                        min.node.size = 5L, seed = 1L)
    pred <- stats::predict(m, data = df_te)$predictions
  } else {
    m <- stats::lm(y ~ ., data = df_tr)
    pred <- stats::predict(m, newdata = df_te)
  }
  rmse <- sqrt(mean((pred - y_test)^2))
  r2 <- 1 - mean((pred - y_test)^2) / mean((y_test - mean(y_train))^2)
  c(rmse = rmse, r2 = r2)
}

landmark_matrix <- function(dat, biomarker_columns) {
  v <- dat$visits
  lm_rows <- do.call(rbind, lapply(split(v, v$id), function(sv) {
    sv[which.min(abs(sv$lag)), , drop = FALSE]
  }))
  lm_rows <- lm_rows[match(dat$subjects$id, lm_rows$id), , drop = FALSE]
  as.matrix(lm_rows[, biomarker_columns, drop = FALSE])
}

full_history_matrix <- function(dat, biomarker_columns) {
  v <- dat$visits
  parts <- lapply(sort(unique(v$visit)), function(pos) {
    sv <- v[v$visit == pos, , drop = FALSE]
    sv <- sv[match(dat$subjects$id, sv$id), , drop = FALSE]
    z <- as.matrix(sv[, biomarker_columns, drop = FALSE])
    colnames(z) <- paste0("v", pos, "_", biomarker_columns)
    z
  })
  z <- do.call(cbind, parts)
  z[!is.finite(z)] <- 0
  z
}

mcse <- function(z) if (length(z) < 2L) NA_real_ else stats::sd(z) / sqrt(length(z))

rows <- lapply(SCENARIOS_ALL, function(scn) {
  sc <- get_scenario(scn)
  reps <- lapply(seq_len(N_REP), function(r) {
    s <- SEED + 1000L * r
    train <- generate_dataset_fixed(N, scn, seed = s)
    test  <- generate_dataset_fixed(N, scn, seed = s + 1L)
    all_markers <- bio_columns(train$visits)
    b1_lm <- fit_predict_rmse(
      landmark_matrix(train, "B1"), train$subjects$A_star,
      landmark_matrix(test, "B1"), test$subjects$A_star
    )
    b1_fh <- fit_predict_rmse(
      full_history_matrix(train, "B1"), train$subjects$A_star,
      full_history_matrix(test, "B1"), test$subjects$A_star
    )
    all_lm <- fit_predict_rmse(
      landmark_matrix(train, all_markers), train$subjects$A_star,
      landmark_matrix(test, all_markers), test$subjects$A_star
    )
    all_fh <- fit_predict_rmse(
      full_history_matrix(train, all_markers), train$subjects$A_star,
      full_history_matrix(test, all_markers), test$subjects$A_star
    )
    c(
      rmse_b1_lm = unname(b1_lm["rmse"]), r2_b1_lm = unname(b1_lm["r2"]),
      rmse_b1_fh = unname(b1_fh["rmse"]), r2_b1_fh = unname(b1_fh["r2"]),
      rmse_all_lm = unname(all_lm["rmse"]), r2_all_lm = unname(all_lm["r2"]),
      rmse_all_fh = unname(all_fh["rmse"]), r2_all_fh = unname(all_fh["r2"])
    )
  })
  M <- do.call(rbind, reps)

  data.frame(
    scenario = scn,
    sigma_eps = sc$sigma_eps,
    signal_amp = sc$signal_amp,
    sigma_bio = sc$sigma_bio,
    n_biomarkers = sc$n_biomarkers,
    schedule = scenario_schedule(scn),
    n_rep = N_REP,
    rmse_b1_landmark = mean(M[, "rmse_b1_lm"]),
    rmse_b1_landmark_mcse = mcse(M[, "rmse_b1_lm"]),
    r2_b1_landmark = mean(M[, "r2_b1_lm"]),
    r2_b1_landmark_mcse = mcse(M[, "r2_b1_lm"]),
    rmse_b1_full_history = mean(M[, "rmse_b1_fh"]),
    rmse_b1_full_history_mcse = mcse(M[, "rmse_b1_fh"]),
    r2_b1_full_history = mean(M[, "r2_b1_fh"]),
    r2_b1_full_history_mcse = mcse(M[, "r2_b1_fh"]),
    rmse_all_landmark = mean(M[, "rmse_all_lm"]),
    rmse_all_landmark_mcse = mcse(M[, "rmse_all_lm"]),
    r2_all_landmark = mean(M[, "r2_all_lm"]),
    r2_all_landmark_mcse = mcse(M[, "r2_all_lm"]),
    rmse_all_full_history = mean(M[, "rmse_all_fh"]),
    rmse_all_full_history_mcse = mcse(M[, "rmse_all_fh"]),
    r2_all_full_history = mean(M[, "r2_all_fh"]),
    r2_all_full_history_mcse = mcse(M[, "r2_all_fh"]),
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)
numeric_cols <- vapply(out, is.numeric, logical(1))
out[, numeric_cols] <- round(out[, numeric_cols], 4)
utils::write.csv(out, OUT, row.names = FALSE)
print(out[, c(
  "scenario", "sigma_eps", "rmse_b1_landmark", "rmse_b1_full_history",
  "rmse_all_landmark", "rmse_all_full_history"
)], row.names = FALSE)
message("\nInterpretation: B1 results are relevant to the current comparator ",
        "feature set;\nall-marker results measure information available to SILK ",
        "and must not be\nattributed to single-marker comparators. Small RMSE ",
        "shows achievability for\nthe diagnostic random forest, not an ",
        "information-theoretic limit or proof that\na fitted survival comparator ",
        "actually recovers the latent stage.")
message("Wrote: ", normalizePath(OUT, mustWork = FALSE))
