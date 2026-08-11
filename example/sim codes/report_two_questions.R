# =============================================================================
# report_two_questions.R
#
# Reporting for the SILK simulation, structured as the two empirical questions
# the manuscript poses. Neither replaces the other.
#
#   Q1  Does registration correct the clock?
#       Hold the feature map g fixed and change ONLY the age coordinate:
#         recorded  ->  SILK-calibrated  ->  latent oracle.
#       Primary map: residual-time Cox on g(age coordinate, invariant history)
#       -- the manuscript's estimator. The attained-age triple is reported
#       separately and labelled a supplementary clock-isolation analysis.
#
#   Q2  Does the complete SILK procedure beat procedures a practitioner could
#       actually run? Prespecified competitor pool, reported honestly.
#
# Estimand note. SILK's consistent estimand conditions on the CALIBRATED age and
# the shift-invariant history, not on the latent age. Primary performance is
# assessed by held-out proper scores (IBS, calibration). Shift and stage RMSE
# against the true shift appear only in Section S as latent-resolution
# diagnostics. Profile separation, boundary frequency and multistart spread
# assess identifiability and optimization; they do not replace predictive
# validation.
#
# Usage:
#   Rscript report_two_questions.R <results_dir> [output_dir] [--rebuild-cache]
#
# The first run consolidates the per-task CSV files into one replication-level
# file under <output_dir>/cache/ and every later run reads only that file.
# =============================================================================

suppressPackageStartupMessages({
  library(stats)
  library(utils)
})

# -------------------------------------------------------------- configuration

args <- commandArgs(trailingOnly = TRUE)
flags <- grep("^--", args, value = TRUE)
positional <- setdiff(args, flags)

RESULTS_DIR <- if (length(positional) >= 1L) positional[[1]] else "."
OUT_DIR <- if (length(positional) >= 2L) positional[[2]] else
  file.path(RESULTS_DIR, "report_two_questions")
REBUILD_CACHE <- "--rebuild-cache" %in% flags
# Completeness is an error by default. --allow-incomplete downgrades it to a
# warning and is for exploratory use only; it must never produce a paper table.
ALLOW_INCOMPLETE <- "--allow-incomplete" %in% flags

CACHE_DIR <- file.path(OUT_DIR, "cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

PRIMARY_N_TRAIN <- as.integer(Sys.getenv("REPORT_N_TRAIN", unset = "400"))
EXPECTED_N_REP <- as.integer(Sys.getenv("REPORT_N_REP", unset = "50"))
N_BOOT <- as.integer(Sys.getenv("REPORT_N_BOOT", unset = "10000"))
BOOT_SEED <- as.integer(Sys.getenv("REPORT_BOOT_SEED", unset = "20260530"))

# ---------------------------------------------------------------------------
# PRESPECIFIED analysis roles. These must match cfg.R. They are written out here
# rather than inferred from whatever methods happen to be present, so that the
# leaderboard cannot silently absorb an ablation, an oracle bound, or a sibling
# SILK kernel, and so that a missing primary method is an error rather than a
# silent fallback.
# ---------------------------------------------------------------------------

# Deliberate, logged override for analysing the legacy 8.9.2026 freeze, which
# has no SILK-History. There is no SILENT fallback: an absent primary method is
# an error unless the analyst names the substitute explicitly here.
PRIMARY_SILK_METHOD <- Sys.getenv("REPORT_PRIMARY_SILK", unset = "SILK-History")
if (!identical(PRIMARY_SILK_METHOD, "SILK-History")) {
  message("NOTE: primary SILK method overridden to '", PRIMARY_SILK_METHOD,
          "' via REPORT_PRIMARY_SILK. This is not the manuscript's estimator; ",
          "use only for legacy-freeze analysis.")
}

TRIPLES <- list(
  primary = list(
    label = "PRIMARY: residual-time Cox, g = (age coordinate, invariant history)",
    recorded = "Cox-History-Recorded",
    silk     = PRIMARY_SILK_METHOD,
    oracle   = "Oracle-History-Latent-Age",
    required = TRUE
  ),
  clock_isolation = list(
    label = "SUPPLEMENTARY: attained-age Cox, g = (X1, X2)  [clock isolation]",
    recorded = "Cox-SameFeature-Recorded",
    silk     = "SILK",
    oracle   = "Oracle-Latent-Age",
    required = FALSE
  )
)

# Complete procedures a practitioner could run. Fixed in advance.
PRIMARY_COMPETITORS <- c(
  "Cox-History-Recorded",
  "Landmark-Recorded",
  "MMLM-Recorded",
  "JM-Recorded",
  "RSF-Observed",
  "DeepSurv-Observed",
  "TimeError-Integrated-Landmark"
)

# Same method under a different kernel or survival layer; never a "competitor".
SILK_FAMILY_METHODS <- c(
  "SILK", "SILK-Laplace", "SILK-Matern32",
  "SILK-History", "SILK-History-Laplace", "SILK-History-Matern32",
  "Beran-SILK"
)

# Kernel-sensitivity arm. Same survival layer and feature map by construction,
# so this contrast varies the kernel alone.
KERNEL_ARM <- c("SILK-History", "SILK-History-Laplace", "SILK-History-Matern32")

# ------------------------------------------------------------- consolidation

read_family_fast <- function(dir, family) {
  files <- list.files(dir, pattern = paste0("^", family, "_task_[0-9]+\\.csv$"),
                      full.names = TRUE)
  if (!length(files)) return(NULL)
  message("  ", family, ": ", length(files), " files in ", basename(dir))
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::rbindlist(
      lapply(files, function(f) {
        d <- data.table::fread(f, showProgress = FALSE)
        d[, task_id := as.integer(sub(".*_task_([0-9]+)\\.csv$", "\\1", f))]
        d
      }),
      use.names = TRUE, fill = TRUE
    )
    return(as.data.frame(dt))
  }
  parts <- lapply(files, function(f) {
    d <- try(read.csv(f, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(d, "try-error") || !nrow(d)) return(NULL)
    d$task_id <- as.integer(sub(".*_task_([0-9]+)\\.csv$", "\\1", f))
    d
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(NULL)
  do.call(rbind, parts)
}

build_cache <- function(family, cache_file) {
  message("Consolidating ", family, " (one-time) ...")
  base <- read_family_fast(file.path(RESULTS_DIR, "raw"), family)
  if (is.null(base)) stop("No ", family, " files under raw/.", call. = FALSE)
  repair_dir <- file.path(RESULTS_DIR, "raw_")
  if (dir.exists(repair_dir)) {
    repaired <- read_family_fast(repair_dir, family)
    if (!is.null(repaired) && nrow(repaired)) {
      common <- intersect(names(base), names(repaired))
      k <- function(d) paste(d$task_id, d$method, sep = "|")
      base <- base[!(k(base) %in% k(repaired)), common, drop = FALSE]
      base <- rbind(base, repaired[, common, drop = FALSE])
      message("  overlaid ", nrow(repaired), " repaired rows")
    }
  }
  saveRDS(base, cache_file)
  base
}

load_family <- function(family) {
  cache_file <- file.path(CACHE_DIR, paste0(family, ".rds"))
  if (!REBUILD_CACHE && file.exists(cache_file)) {
    message("Reading cached ", family, " from ", basename(cache_file))
    return(readRDS(cache_file))
  }
  build_cache(family, cache_file)
}

metrics <- load_family("prediction_metrics")

parts <- strsplit(as.character(metrics$time_grid_setting), ":", fixed = TRUE)
metrics$phase <- vapply(parts, `[`, character(1), 1L)
metrics$scenario <- vapply(parts, `[`, character(1), 2L)
metrics$schedule <- vapply(parts, `[`, character(1), 3L)

primary <- metrics[metrics$phase == "primary" &
                     metrics$n_train_setting == PRIMARY_N_TRAIN, , drop = FALSE]
if (!nrow(primary)) {
  stop("No primary rows at n_train = ", PRIMARY_N_TRAIN, ".", call. = FALSE)
}

# ------------------------------------------------------------- validation ---

problems <- character(0)
note <- function(...) problems <<- c(problems, paste0(...))

key <- paste(primary$scenario, primary$replicate_id, primary$method, sep = "|")
dup <- key[duplicated(key)]
if (length(dup)) {
  note("Duplicate scenario-replicate-method keys: ", length(unique(dup)),
       " (e.g. ", paste(utils::head(unique(dup), 3L), collapse = "; "), ")")
}

scenarios <- sort(unique(primary$scenario))
available <- sort(unique(primary$method))

required_methods <- unique(c(
  unlist(lapply(TRIPLES[vapply(TRIPLES, function(z) isTRUE(z$required), logical(1))],
                function(z) c(z$recorded, z$silk, z$oracle))),
  PRIMARY_COMPETITORS
))
missing_required <- setdiff(required_methods, available)
if (length(missing_required)) {
  note("Prespecified methods absent from the results: ",
       paste(missing_required, collapse = ", "),
       ". There is no fallback: rerun the simulation with the correct roster.")
}

for (m in intersect(required_methods, available)) {
  counts <- table(primary$scenario[primary$method == m])
  short <- names(counts)[counts != EXPECTED_N_REP]
  absent <- setdiff(scenarios, names(counts))
  if (length(short)) {
    note("Method ", m, " does not have exactly ", EXPECTED_N_REP,
         " replicates in: ", paste(short, collapse = ", "))
  }
  if (length(absent)) {
    note("Method ", m, " missing entirely in: ", paste(absent, collapse = ", "))
  }
}

if (length(problems)) {
  msg <- paste0("Completeness / integrity failures:\n  - ",
                paste(problems, collapse = "\n  - "))
  if (ALLOW_INCOMPLETE) {
    warning(msg, "\nProceeding because --allow-incomplete was set. ",
            "Do NOT use this output in the manuscript.", call. = FALSE)
  } else {
    stop(msg,
         "\n\nFix the run, or pass --allow-incomplete for exploratory use only.",
         call. = FALSE)
  }
} else {
  message("Completeness check passed: ", length(scenarios), " scenarios x ",
          EXPECTED_N_REP, " replicates x ", length(required_methods),
          " prespecified methods.")
}

wide <- reshape(
  primary[, c("scenario", "replicate_id", "method", "integrated_brier_score")],
  idvar = c("scenario", "replicate_id"), timevar = "method", direction = "wide"
)
names(wide) <- sub("^integrated_brier_score\\.", "", names(wide))
wide <- wide[order(wide$scenario, wide$replicate_id), , drop = FALSE]

# --------------------------------------------------------------- bootstrap ---

set.seed(BOOT_SEED)

boot_index <- function(n, B = N_BOOT) {
  matrix(sample.int(n, n * B, replace = TRUE), nrow = B)
}

boot_mean_diff <- function(a, b) {
  keep <- is.finite(a) & is.finite(b)
  a <- a[keep]
  b <- b[keep]
  if (!length(a)) return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  idx <- boot_index(length(a))
  draws <- rowMeans(matrix(a[idx], nrow = nrow(idx)) -
                      matrix(b[idx], nrow = nrow(idx)))
  c(estimate = mean(a - b),
    lower = unname(quantile(draws, 0.025)),
    upper = unname(quantile(draws, 0.975)))
}

# ------------------------------------------------------------------- Q1 ------

q1_table <- function(map) {
  need <- c(map$recorded, map$silk, map$oracle)
  if (!all(need %in% names(wide))) return(NULL)
  rows <- lapply(scenarios, function(s) {
    g <- wide[wide$scenario == s, , drop = FALSE]
    g <- g[complete.cases(g[, need]), , drop = FALSE]
    if (nrow(g) != EXPECTED_N_REP && !ALLOW_INCOMPLETE) {
      stop("Scenario ", s, " has ", nrow(g), " complete replicates for the ",
           "triple; expected ", EXPECTED_N_REP, ".", call. = FALSE)
    }
    rec <- g[[map$recorded]]; sil <- g[[map$silk]]; orc <- g[[map$oracle]]

    # The oracle gap is the denominator of "percent closed". When recorded and
    # oracle predictions coincide the ratio is undefined; when they are close
    # the ratio is unstable. `denominator_separated` means "far enough from zero
    # for the ratio to be reported", NOT "mathematically defined".
    recorded_oracle <- boot_mean_diff(rec, orc)
    separated <- !(recorded_oracle[["lower"]] <= 0 &&
                     recorded_oracle[["upper"]] >= 0)

    recorded_silk <- boot_mean_diff(rec, sil)
    silk_oracle <- boot_mean_diff(sil, orc)
    idx <- boot_index(nrow(g))
    ratio <- NULL
    if (separated) {
      nb <- rowMeans(matrix(rec[idx], nrow = nrow(idx)) -
                       matrix(sil[idx], nrow = nrow(idx)))
      db <- rowMeans(matrix(rec[idx], nrow = nrow(idx)) -
                       matrix(orc[idx], nrow = nrow(idx)))
      r <- ifelse(abs(db) > 1e-12, 100 * nb / db, NA_real_)
      ratio <- c(100 * mean(rec - sil) / mean(rec - orc),
                 unname(quantile(r, 0.025, na.rm = TRUE)),
                 unname(quantile(r, 0.975, na.rm = TRUE)))
    }

    data.frame(
      scenario = s, n_rep = nrow(g),
      ibs_recorded = mean(rec), ibs_silk = mean(sil), ibs_oracle = mean(orc),
      silk_minus_oracle_x1e3 = 1000 * silk_oracle[["estimate"]],
      silk_minus_oracle_lo_x1e3 = 1000 * silk_oracle[["lower"]],
      silk_minus_oracle_hi_x1e3 = 1000 * silk_oracle[["upper"]],
      recorded_minus_silk_x1e3 = 1000 * recorded_silk[["estimate"]],
      recorded_minus_silk_lo_x1e3 = 1000 * recorded_silk[["lower"]],
      recorded_minus_silk_hi_x1e3 = 1000 * recorded_silk[["upper"]],
      recorded_minus_oracle_x1e3 = 1000 * recorded_oracle[["estimate"]],
      denominator_separated = separated,
      pct_gap_closed = if (separated) ratio[[1]] else NA_real_,
      pct_lo = if (separated) ratio[[2]] else NA_real_,
      pct_hi = if (separated) ratio[[3]] else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

for (nm in names(TRIPLES)) {
  map <- TRIPLES[[nm]]
  tab <- q1_table(map)
  if (is.null(tab)) {
    if (isTRUE(map$required)) {
      stop("Required triple '", nm, "' is absent from the results.", call. = FALSE)
    }
    message("Q1 [", nm, "]: not present; skipped.")
    next
  }
  write.csv(tab, file.path(OUT_DIR, paste0("q1_", nm, ".csv")), row.names = FALSE)
  message("\n== Q1  ", map$label, " ==")
  show <- data.frame(
    scenario = tab$scenario,
    silk_minus_oracle_x1e3 = round(tab$silk_minus_oracle_x1e3, 3),
    paired_ci = sprintf(
      "[%.3f, %.3f]",
      tab$silk_minus_oracle_lo_x1e3,
      tab$silk_minus_oracle_hi_x1e3
    ),
    recorded_minus_silk_x1e3 = round(tab$recorded_minus_silk_x1e3, 3),
    stringsAsFactors = FALSE
  )
  print(show, row.names = FALSE)
}
message(
  "\nQ1 primary quantity is paired IBS(SILK) - IBS(Oracle); values close to zero",
  "\nwith their paired 95% confidence intervals directly answer the oracle-gap question.",
  "\nRecorded-minus-SILK is positive when SILK improves on plain recorded time."
)

# ------------------------------------------------------------------- Q2 ------

competitors <- intersect(PRIMARY_COMPETITORS, names(wide))
if (!PRIMARY_SILK_METHOD %in% names(wide)) {
  stop("Prespecified primary method '", PRIMARY_SILK_METHOD,
       "' is absent. No fallback is permitted.", call. = FALSE)
}

q2_rows <- list()
verdict_rows <- list()
for (s in scenarios) {
  g <- wide[wide$scenario == s, , drop = FALSE]
  ref <- g[[PRIMARY_SILK_METHOD]]
  for (m in competitors) {
    d <- boot_mean_diff(g[[m]], ref)
    status <- if (d[["lower"]] > 0) {
      "SILK win"
    } else if (d[["upper"]] < 0) {
      "SILK loss"
    } else {
      "statistical tie"
    }
    q2_rows[[length(q2_rows) + 1L]] <- data.frame(
      scenario = s, method = m, mean_ibs = mean(g[[m]]),
      delta_vs_silk_x1e3 = 1000 * d[["estimate"]],
      lo_x1e3 = 1000 * d[["lower"]], hi_x1e3 = 1000 * d[["upper"]],
      result = status,
      stringsAsFactors = FALSE
    )
  }

  # Selection-aware inference. The "best feasible competitor" is chosen from the
  # data, so an ordinary paired interval for the winner understates uncertainty.
  # The best competitor is therefore RESELECTED inside every bootstrap draw.
  complete <- complete.cases(g[, c(PRIMARY_SILK_METHOD, competitors), drop = FALSE])
  g_best <- g[complete, , drop = FALSE]
  if (nrow(g_best) != EXPECTED_N_REP && !ALLOW_INCOMPLETE) {
    stop("Scenario ", s, " has ", nrow(g_best),
         " complete replicates for Q2; expected ", EXPECTED_N_REP, ".",
         call. = FALSE)
  }
  ref_best <- g_best[[PRIMARY_SILK_METHOD]]
  idx <- boot_index(nrow(g_best))
  M <- vapply(competitors, function(m) g_best[[m]], numeric(nrow(g_best)))
  draws <- numeric(nrow(idx))
  picked <- character(nrow(idx))
  for (b in seq_len(nrow(idx))) {
    i <- idx[b, ]
    means <- colMeans(M[i, , drop = FALSE])
    w <- which.min(means)
    picked[b] <- competitors[w]
    draws[b] <- mean(M[i, w] - ref_best[i])
  }
  observed_best <- competitors[which.min(colMeans(M))]
  lo <- unname(quantile(draws, 0.025)); hi <- unname(quantile(draws, 0.975))
  verdict_rows[[length(verdict_rows) + 1L]] <- data.frame(
    scenario = s,
    best_feasible_competitor = observed_best,
    selection_stability = round(mean(picked == observed_best), 3),
    delta_x1e3 = round(1000 * mean(M[, observed_best] - ref_best), 3),
    ci_selection_aware = sprintf("[%.3f, %.3f]", 1000 * lo, 1000 * hi),
    verdict = if (lo > 0) "SILK win" else if (hi < 0) "SILK loss" else "statistical tie",
    stringsAsFactors = FALSE
  )
}
q2 <- do.call(rbind, q2_rows)
verdict <- do.call(rbind, verdict_rows)
q2_summary <- do.call(rbind, lapply(split(q2, q2$scenario), function(z) {
  losses <- z$method[z$result == "SILK loss"]
  data.frame(
    scenario = z$scenario[1L],
    n_comparators = nrow(z),
    silk_wins = sum(z$result == "SILK win"),
    statistical_ties = sum(z$result == "statistical tie"),
    silk_losses = length(losses),
    win_or_tie_against_all = length(losses) == 0L,
    losses_to = if (length(losses)) paste(losses, collapse = "; ") else "",
    stringsAsFactors = FALSE
  )
}))
write.csv(q2, file.path(OUT_DIR, "q2_feasible_leaderboard.csv"), row.names = FALSE)
write.csv(verdict, file.path(OUT_DIR, "q2_verdict.csv"), row.names = FALSE)
write.csv(q2_summary, file.path(OUT_DIR, "q2_win_tie_summary.csv"), row.names = FALSE)
message("\n== Q2  complete feasible procedures (reference: ", PRIMARY_SILK_METHOD, ") ==")
print(q2_summary, row.names = FALSE)
message("Intervals reselect the best competitor within each bootstrap draw, so",
        "\nthey account for choosing the winner from the data. selection_stability",
        "\nis the fraction of draws in which the observed winner won again.")

# ------------------------------------------------------- kernel sensitivity --

kernel_arm <- intersect(KERNEL_ARM, names(wide))
if (length(kernel_arm) >= 2L) {
  ref <- kernel_arm[[1]]
  krows <- lapply(scenarios, function(s) {
    g <- wide[wide$scenario == s, , drop = FALSE]
    do.call(rbind, lapply(setdiff(kernel_arm, ref), function(m) {
      d <- boot_mean_diff(g[[m]], g[[ref]])
      data.frame(scenario = s, kernel = m, reference = ref,
                 delta_x1e3 = 1000 * d[["estimate"]],
                 lo_x1e3 = 1000 * d[["lower"]], hi_x1e3 = 1000 * d[["upper"]],
                 stringsAsFactors = FALSE)
    }))
  })
  write.csv(do.call(rbind, krows),
            file.path(OUT_DIR, "kernel_sensitivity.csv"), row.names = FALSE)
  message("\nKernel sensitivity written (identical survival layer, feature map, ",
          "folds and seeds; kernel is the only difference).")
} else {
  message("\nKernel sensitivity skipped: fewer than two kernel variants present.")
}

# ---------------------------------------- S: supplementary diagnostics only --

reg <- try(load_family("registration_diagnostics"), silent = TRUE)
if (!inherits(reg, "try-error") && !is.null(reg)) {
  r <- reg[reg$phase == "primary" & reg$n_train == PRIMARY_N_TRAIN &
             reg$sample_role == "test_heldout", , drop = FALSE]
  wanted <- c("stage_rmse_improvement", "shift_rmse", "boundary_fraction",
              "median_gap_q1", "median_gap_q2", "clock_signal_r2",
              "longitudinal_kernel_signal", "n_successful_starts",
              "objective_spread", "second_best_gap",
              "fold_mean_objective_spread")
  wanted <- intersect(wanted, names(r))
  agg <- aggregate(r[, wanted, drop = FALSE],
                   by = list(scenario = r$scenario,
                             biomarker_kernel = r$biomarker_kernel),
                   FUN = function(z) mean(z, na.rm = TRUE))
  write.csv(agg, file.path(OUT_DIR, "s_supplementary_diagnostics.csv"),
            row.names = FALSE)
  message("\n== S  supplementary diagnostics ==")
  message("Identifiability / optimization (estimand-relevant, but NOT predictive ",
          "validation):\n  median_gap_q1, median_gap_q2, boundary_fraction, ",
          "n_successful_starts,\n  objective_spread, second_best_gap.")
  message("Latent-resolution only (NOT the estimand; debugging and the appendix ",
          "resolution\n  regime): shift_rmse, stage_rmse_improvement.")
}

message("\nWrote report to: ", normalizePath(OUT_DIR, mustWork = FALSE))
