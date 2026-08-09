# =============================================================================
# data-generation.R
# Data generating mechanism (from dgm.R)
# =============================================================================

#' @keywords internal
logistic <- function(x) 1 / (1 + exp(-x))

#' @keywords internal
r_common_shift <- function(n, sigma_eps, eps_mean = 0, eps_type = "normal") {
  if (sigma_eps <= 0) return(rep(eps_mean, n))
  if (eps_type == "normal") {
    eps_mean + stats::rnorm(n, sd = sigma_eps)
  } else if (eps_type == "mixture") {
    tail_prob <- 0.20
    core_sd <- 0.35
    tail_sd <- 6.00
    is_tail <- stats::rbinom(n, size = 1L, prob = tail_prob)
    z <- stats::rnorm(n, sd = ifelse(is_tail == 1L, tail_sd, core_sd))
    mix_sd <- sqrt((1 - tail_prob) * core_sd^2 + tail_prob * tail_sd^2)
    eps_mean + z / mix_sd * sigma_eps
  } else if (eps_type == "laplace") {
    z <- (stats::rexp(n) - stats::rexp(n)) / sqrt(2)
    eps_mean + z * sigma_eps
  } else if (eps_type == "t3") {
    df <- 2.5
    eps_mean + stats::rt(n, df = df) / sqrt(df / (df - 2)) * sigma_eps
  } else if (eps_type == "asymmetric") {
    z <- stats::rchisq(n, df = 1) - 1
    eps_mean + z / stats::sd(z) * sigma_eps
  } else {
    stop("Unknown eps_type: ", eps_type, call. = FALSE)
  }
}

#' @keywords internal
true_survival_curve <- function(u, eta_risk) {
  SURV_LAMBDA_D <- silk_opt("SURV_LAMBDA_D")
  SURV_KAPPA_D <- silk_opt("SURV_KAPPA_D")
  exp(-((u / SURV_LAMBDA_D) ^ SURV_KAPPA_D) * exp(eta_risk))
}

#' @keywords internal
rweibull_ph <- function(eta_risk) {
  SURV_LAMBDA_D <- silk_opt("SURV_LAMBDA_D")
  SURV_KAPPA_D <- silk_opt("SURV_KAPPA_D")
  u <- stats::runif(length(eta_risk))
  SURV_LAMBDA_D * (-log1p(-u) * exp(-eta_risk)) ^ (1 / SURV_KAPPA_D)
}

#' @keywords internal
sc_val <- function(sc, nm, default) {
  z <- sc[[nm]][1]
  if (length(z) == 0L || is.na(z)) default else z
}

#' @keywords internal
marker_mean_matrix <- function(age, X1, X2, U, p, sc) {
  A_STAR_CENTER <- silk_opt("A_STAR_CENTER")
  DEFAULT_SIGNAL_AMP <- silk_opt("DEFAULT_SIGNAL_AMP")
  DEFAULT_U_BIO_COEF <- silk_opt("DEFAULT_U_BIO_COEF")

  n <- length(age)
  z <- (age - A_STAR_CENTER) / 2.5
  M <- matrix(0, nrow = n, ncol = p)

  signal <- sc$biomarker_signal[1]
  amp <- sc_val(sc, "signal_amp", if (signal == "weak") 0.12 else DEFAULT_SIGNAL_AMP)
  ucoef <- sc_val(sc, "u_bio_coef", DEFAULT_U_BIO_COEF)

  base_list <- list(
    2.2 * tanh(z),
    -1.8 * tanh((age - 29) / 2.0),
    1.5 * sin((age - 24) * pi / 8),
    1.2 * cos((age - 25) * pi / 7),
    0.35 * (age - A_STAR_CENTER)^2 - 2.0
  )

  for (k in seq_len(p)) {
    bk <- base_list[[((k - 1L) %% length(base_list)) + 1L]]
    M[, k] <- amp * bk + 0.18 * X1 + 0.10 * X2 + ucoef * U
  }
  M
}

#' @keywords internal
r_biomarkers <- function(age, X1, X2, U, sc, p) {
  A_STAR_CENTER <- silk_opt("A_STAR_CENTER")
  A_STAR_SD <- silk_opt("A_STAR_SD")
  DIST_SD_BASE <- silk_opt("DIST_SD_BASE")
  DIST_SD_SLOPE <- silk_opt("DIST_SD_SLOPE")
  SIGMA_BIO_MEAN <- silk_opt("SIGMA_BIO_MEAN")
  SIGMA_BIO_WEAK <- silk_opt("SIGMA_BIO_WEAK")
  DEFAULT_U_BIO_COEF <- silk_opt("DEFAULT_U_BIO_COEF")

  n <- length(age)
  signal <- sc$biomarker_signal[1]

  if (signal %in% c("distribution", "distribution_strong")) {
    ucoef <- sc_val(sc, "u_bio_coef", 0.10)
    mean_mat <- matrix(0.08 * X1 + 0.05 * X2 + ucoef * U, nrow = n, ncol = p)
    stage_z <- (age - A_STAR_CENTER) / A_STAR_SD
    sd_age <- sc_val(sc, "dist_sd_base", DIST_SD_BASE) *
      exp(sc_val(sc, "dist_sd_slope", DIST_SD_SLOPE) * stage_z / 2)
    sd_age <- pmax(sd_age, 0.10)

    Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)

    tail_prob <- logistic(-1.35 + 1.15 * stage_z)
    tail_scale <- ifelse(stats::rbinom(n, size = 1L, prob = tail_prob) == 1L, 2.75, 1.0)
    Z <- Z * matrix(tail_scale, nrow = n, ncol = p, byrow = FALSE)

    common_sd <- 0.10 + 0.22 * logistic(stage_z)
    common <- stats::rnorm(n) * common_sd
    out <- mean_mat + Z * matrix(sd_age, nrow = n, ncol = p, byrow = FALSE)
    if (p >= 4L) {
      out[, seq(1, p, by = 4)] <- out[, seq(1, p, by = 4), drop = FALSE] + common
    }
  } else if (signal == "weak") {
    mean_mat <- marker_mean_matrix(age, X1, X2, U, p, sc)
    out <- mean_mat + matrix(stats::rnorm(n * p, sd = sc_val(sc, "sigma_bio", SIGMA_BIO_WEAK)), nrow = n, ncol = p)
  } else {
    mean_mat <- marker_mean_matrix(age, X1, X2, U, p, sc)
    out <- mean_mat + matrix(stats::rnorm(n * p, sd = sc_val(sc, "sigma_bio", SIGMA_BIO_MEAN)), nrow = n, ncol = p)
  }

  colnames(out) <- paste0("B", seq_len(p))
  out
}

#' @keywords internal
make_visits <- function(subjects, lags, sc) {
  n <- nrow(subjects)
  rows <- vector("list", n)
  p <- as.integer(sc$n_biomarkers)

  for (i in seq_len(n)) {
    lags_i <- lags
    # Preserve the nominal schedule position before irregular timing and
    # missingness are applied. Renumbering retained visits would incorrectly
    # align different biological positions in the registration template.
    visit_position <- seq_along(lags_i)
    if (isTRUE(sc$irregular)) {
      jitter <- stats::runif(length(lags_i), min = -0.25, max = 0.25)
      jitter[lags_i == 0] <- 0
      lags_i <- pmax(lags_i + jitter, 0)
      lags_i[lags == 0] <- 0
    }

    keep <- rep(TRUE, length(lags_i))
    if (sc$missing_rate > 0 && length(lags_i) > 1L) {
      pre <- which(lags_i > 0)
      keep[pre] <- stats::runif(length(pre)) > sc$missing_rate
      keep[lags_i == 0] <- TRUE
      if (sum(keep) < 2L && length(pre) > 0L) keep[sample(pre, 1)] <- TRUE
    }

    lags_i <- lags_i[keep]
    visit_position <- visit_position[keep]
    visit_age_star <- subjects$A_star[i] - lags_i
    visit_age_obs <- visit_age_star + subjects$eps[i]
    B <- r_biomarkers(
      age = visit_age_star,
      X1 = rep(subjects$X1[i], length(lags_i)),
      X2 = rep(subjects$X2[i], length(lags_i)),
      U = rep(subjects$U_latent[i], length(lags_i)),
      sc = sc,
      p = p
    )

    rows[[i]] <- cbind(
      data.frame(
        id = subjects$id[i],
        visit = visit_position,
        lag = lags_i,
        A_star_il = visit_age_star,
        A_obs_il = visit_age_obs,
        X1 = subjects$X1[i],
        X2 = subjects$X2[i],
        stringsAsFactors = FALSE
      ),
      B
    )
  }

  do.call(rbind, rows)
}

#' Generate a simulated dataset
#'
#' Generates training or test data under a specified error scenario and visit
#' schedule for the SILK simulation framework.
#'
#' @param n Integer. Number of subjects.
#' @param scenario_name Character string. Name of the error scenario.
#' @param schedule_name Character string. Visit schedule name (default: scenario default).
#' @param seed Integer or NULL. Random seed for reproducibility.
#' @return A list with elements:
#'   \describe{
#'     \item{subjects}{Data frame of subject-level information.}
#'     \item{visits}{Data frame of visit-level biomarker measurements.}
#'     \item{scenario}{The scenario name used.}
#'     \item{schedule}{The schedule name used.}
#'     \item{scenario_row}{The full scenario definition row.}
#'   }
#' @export
#' @examples
#' dat <- generate_dataset_fixed(100, "mean_moderate", seed = 42)
#' head(dat$subjects)
#' head(dat$visits)
generate_dataset_fixed <- function(n, scenario_name, schedule_name = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sc <- get_scenario(scenario_name)
  if (is.null(schedule_name) || is.na(schedule_name) || !nzchar(schedule_name)) schedule_name <- scenario_schedule(scenario_name)
  VISIT_SCHEDULES <- silk_opt("VISIT_SCHEDULES")
  lags <- VISIT_SCHEDULES[[schedule_name]]
  if (is.null(lags)) stop("Unknown schedule: ", schedule_name, call. = FALSE)

  LATENT_AGE_MIN <- silk_opt("LATENT_AGE_MIN")
  LATENT_AGE_MAX <- silk_opt("LATENT_AGE_MAX")
  A_STAR_CENTER <- silk_opt("A_STAR_CENTER")
  A_STAR_SD <- silk_opt("A_STAR_SD")
  SURV_BETA_A <- silk_opt("SURV_BETA_A")
  SURV_BETA_X1 <- silk_opt("SURV_BETA_X1")
  SURV_BETA_X2 <- silk_opt("SURV_BETA_X2")
  SURV_BETA_U <- silk_opt("SURV_BETA_U")
  CENS_RATE <- silk_opt("CENS_RATE")

  X1 <- stats::rnorm(n)
  X2 <- stats::rbinom(n, 1, 0.45)
  U_latent <- stats::rnorm(n)
  A_star <- stats::runif(n, LATENT_AGE_MIN, LATENT_AGE_MAX)
  eps <- r_common_shift(n, sc$sigma_eps, sc$eps_mean, sc$eps_type)
  A_obs <- A_star + eps

  beta_A <- sc_val(sc, "risk_beta_A", SURV_BETA_A)
  beta_U <- sc_val(sc, "risk_beta_U", SURV_BETA_U)
  eta_risk <- beta_A * ((A_star - A_STAR_CENTER) / A_STAR_SD) +
    SURV_BETA_X1 * X1 + SURV_BETA_X2 * X2 + beta_U * U_latent

  D_star <- rweibull_ph(eta_risk)
  G_star <- stats::rexp(n, rate = CENS_RATE)
  U <- pmin(D_star, G_star)
  delta <- as.integer(D_star <= G_star)
  T_star <- A_star + D_star
  C_star <- A_star + G_star
  T_obs <- T_star + eps
  C_obs <- C_star + eps

  subjects <- data.frame(
    id = seq_len(n),
    X1 = X1,
    X2 = X2,
    U_latent = U_latent,
    A_star = A_star,
    eps = eps,
    A_obs = A_obs,
    eta_risk = eta_risk,
    D_star = D_star,
    G_star = G_star,
    T_star = T_star,
    C_star = C_star,
    T_obs = T_obs,
    C_obs = C_obs,
    U = U,
    delta = delta,
    stringsAsFactors = FALSE
  )

  visits <- make_visits(subjects, lags, sc)

  list(
    subjects = subjects,
    visits = visits,
    scenario = scenario_name,
    schedule = schedule_name,
    scenario_row = sc
  )
}

#' Get biomarker column names
#'
#' Returns the names of biomarker columns (B1, B2, ...) in a data frame.
#'
#' @param df Data frame with biomarker columns.
#' @return Character vector of biomarker column names.
#' @export
bio_columns <- function(df) {
  cols <- grep("^B[0-9]+$", names(df), value = TRUE)
  if (length(cols) == 0L) stop("No biomarker columns B1, B2, ... found", call. = FALSE)
  cols
}

#' Build landmark history features
#'
#' Creates summary features from biomarker visit histories for use in survival
#' prediction models.
#'
#' @param subjects Data frame of subjects.
#' @param visits Data frame of visits with biomarker columns.
#' @param max_biomarkers Integer. Maximum number of biomarkers to include as
#'   current values (default 4).
#' @return Data frame with one row per subject containing history features.
#' @export
make_history_features <- function(subjects, visits, max_biomarkers = 4L) {
  ids <- subjects$id
  bcols_all <- bio_columns(visits)
  bcols <- bcols_all[seq_len(min(length(bcols_all), max_biomarkers))]

  out <- data.frame(id = ids, stringsAsFactors = FALSE)
  for (nm in bcols) {
    out[[paste0("current_", nm)]] <- NA_real_
  }
  out$mean_B <- NA_real_
  out$sd_B <- NA_real_
  out$visit_count <- NA_real_
  out$span <- NA_real_
  out$X1 <- subjects$X1[match(ids, subjects$id)]
  out$X2 <- subjects$X2[match(ids, subjects$id)]

  for (k in seq_along(ids)) {
    sv <- visits[visits$id == ids[k], , drop = FALSE]
    sv <- sv[order(sv$lag, decreasing = TRUE), , drop = FALSE]
    lm <- sv[which.min(abs(sv$lag)), , drop = FALSE]
    out$visit_count[k] <- nrow(sv)
    out$span[k] <- max(sv$lag) - min(sv$lag)

    Bmat <- as.matrix(sv[, bcols_all, drop = FALSE])
    out$mean_B[k] <- mean(Bmat)
    out$sd_B[k] <- stats::sd(as.vector(Bmat))

    for (nm in bcols) {
      out[[paste0("current_", nm)]][k] <- lm[[nm]][1]
    }
  }

  out[is.na(out)] <- 0
  out
}

#' @keywords internal
make_true_survival_matrix <- function(subjects, u_grid) {
  mat <- outer(subjects$eta_risk, u_grid, function(eta, u) true_survival_curve(u, eta))
  matrix(mat, nrow = nrow(subjects), ncol = length(u_grid))
}
