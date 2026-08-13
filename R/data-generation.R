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
    # Contaminated normal. The component SDs are specified RELATIVE to each
    # other and the draw is then rescaled to marginal SD sigma_eps, so what
    # matters is the ratio tail_sd / core_sd, not their absolute values.
    # A large ratio concentrates almost all mass at a near-zero shift: with the
    # former (core 0.35, tail 6.00) the rescaled core SD was only
    # 0.35 * sigma_eps / 2.70 = 0.65 years at sigma_eps = 5, leaving 80% of
    # subjects with an essentially correct clock and P(|eps| < 1) = 72%.
    # The current values give a rescaled core SD of 3.10 and tail SD of 9.30 at
    # sigma_eps = 5, i.e. a genuine common shift for every subject with a
    # threefold-worse contaminating component.
    tail_prob <- as.numeric(silk_opt("MIX_TAIL_PROB"))[1L]
    core_sd <- as.numeric(silk_opt("MIX_CORE_SD"))[1L]
    tail_sd <- as.numeric(silk_opt("MIX_TAIL_SD"))[1L]
    is_tail <- stats::rbinom(n, size = 1L, prob = tail_prob)
    z <- stats::rnorm(n, sd = ifelse(is_tail == 1L, tail_sd, core_sd))
    mix_sd <- sqrt((1 - tail_prob) * core_sd^2 + tail_prob * tail_sd^2)
    eps_mean + z / mix_sd * sigma_eps
  } else if (eps_type == "laplace") {
    z <- (stats::rexp(n) - stats::rexp(n)) / sqrt(2)
    eps_mean + z * sigma_eps
  } else if (eps_type == "t3") {
    # Heavy-tailed origin error. Standardising a t to marginal SD sigma_eps
    # inflates the interquartile range and therefore produces a WIDER bulk and
    # only modestly heavier tails than N(0, sigma_eps) -- the opposite of the
    # intended contrast. Matching the interquartile range instead keeps the
    # bulk comparable to N(0, sigma_eps) while leaving genuinely heavy tails:
    # at sigma_eps = 5, P(|eps| > 15) rises from 0.3% under the normal to 5.2%.
    # Note that sigma_eps is then a normal-equivalent scale, not the marginal
    # SD, which is not finite in a useful sense for a t with df = 2.5.
    df <- as.numeric(silk_opt("HEAVY_TAIL_DF"))[1L]
    scale_iqr <- sigma_eps / (stats::qt(0.75, df) / stats::qnorm(0.75))
    eps_mean + stats::rt(n, df = df) * scale_iqr
  } else if (eps_type == "asymmetric") {
    z <- stats::rchisq(n, df = 1) - 1
    eps_mean + z / stats::sd(z) * sigma_eps
  } else {
    stop("Unknown eps_type: ", eps_type, call. = FALSE)
  }
}

#' @keywords internal
attained_age_cumulative_hazard <- function(age, growth = NULL) {
  if (is.null(growth)) growth <- silk_opt("SURV_AGE_GROWTH")
  rate <- as.numeric(silk_opt("SURV_AGE_BASE_RATE"))[1L]
  growth <- as.numeric(growth)[1L]
  center <- as.numeric(silk_opt("A_STAR_CENTER"))[1L]
  if (!is.finite(rate) || rate <= 0 || !is.finite(growth) || growth <= 0) {
    stop("Attained-age Gompertz rate and growth must be positive.", call. = FALSE)
  }
  rate / growth * exp(growth * (as.numeric(age) - center))
}

#' @keywords internal
true_survival_curve <- function(u, start_age, eta_risk, growth = NULL) {
  cumulative_increment <- attained_age_cumulative_hazard(start_age + u, growth) -
    attained_age_cumulative_hazard(start_age, growth)
  exp(-pmax(cumulative_increment, 0) * exp(eta_risk))
}

#' @keywords internal
r_attained_age_ph <- function(start_age, eta_risk, growth = NULL) {
  target <- stats::rexp(length(eta_risk)) * exp(-eta_risk)
  start_hazard <- attained_age_cumulative_hazard(start_age, growth)
  growth <- if (is.null(growth)) silk_opt("SURV_AGE_GROWTH") else growth
  rate <- silk_opt("SURV_AGE_BASE_RATE")
  center <- silk_opt("A_STAR_CENTER")
  stop_age <- center + log((start_hazard + target) * growth / rate) / growth
  pmax(stop_age - start_age, 0)
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

  # Template channels. These must identify a common shift over the WHOLE
  # candidate grid, which is roughly +/-3 sigma_eps and therefore +/-36 at the
  # severe setting. Two earlier channels violated that:
  #
  #   1.5 * sin((age - 24) * pi / 8)   period 16
  #   1.2 * cos((age - 25) * pi / 7)   period 14
  #
  # nearly re-align at a 15-year shift, so with four biomarkers the population
  # discrepancy || g(a) - g(a + d) ||^2 had local minima at d = -15.0 (22% of
  # the maximum loss), -30.4, +14.8 and +30.2 in addition to the truth at 0.
  # This was harmless while the grid stopped at +/-15 -- the competing basin sat
  # exactly on the boundary -- but at +/-36 every subject has four basins.
  # The tanh channels saturate beyond about +/-7 and cannot break the tie.
  #
  # The fifth channel, 0.35 * (age - 28)^2 - 2, was symmetric about the age
  # center and carried SD 9.95 against roughly 1.2 for every other channel, so
  # in the scenarios that include it registration was driven by that one
  # channel alone.
  #
  # Replacements: staggered saturating channels at different centers and widths
  # for g3/g4, and a monotone cubic on a comparable scale for g5. The
  # discrepancy is then strictly increasing in |d| over +/-36 with no spurious
  # minima, and near-zero discrimination improves (normalised loss at |d| = 1
  # goes from 0.022 to 0.014).
  base_list <- list(
    2.2 * tanh(z),
    -1.8 * tanh((age - 29) / 2.0),
    1.5 * tanh((age - 25) / 4.0),
    -1.2 * tanh((age - 31) / 3.0),
    1.4 * ((age - A_STAR_CENTER) / 5) + 0.30 * ((age - A_STAR_CENTER) / 5)^3
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
    # A deliberately weak mean trajectory supplies a stable clock initializer;
    # the principal signal remains the stage-varying, moment-matched shape
    # below, which is what the characteristic-kernel loss exploits.
    mean_mat <- marker_mean_matrix(age, X1, X2, U, p, sc)
    stage_z <- (age - A_STAR_CENTER) / A_STAR_SD
    # Shape-only signal: the stage-dependent innovation moves from Gaussian to
    # symmetric bimodal while preserving mean zero and variance one. Hence raw
    # means and variances do not identify stage, whereas a characteristic
    # kernel mean embedding can distinguish the full distributions.
    shape_weight <- logistic(
      sc_val(sc, "dist_sd_slope", DIST_SD_SLOPE) * stage_z
    )
    use_bimodal <- matrix(
      stats::runif(n * p) < rep(shape_weight, each = p), nrow = n, byrow = TRUE
    )
    separation <- if (signal == "distribution_strong") 0.90 else 0.78
    component <- matrix(sample(c(-1, 1), n * p, replace = TRUE), nrow = n)
    normal <- matrix(stats::rnorm(n * p), nrow = n)
    bimodal <- separation * component + sqrt(1 - separation^2) * normal
    innovation <- ifelse(use_bimodal, bimodal, normal)
    scale <- sc_val(sc, "dist_sd_base", DIST_SD_BASE)
    out <- mean_mat + scale * innovation
  } else if (signal == "stage_null") {
    # Prespecified negative control: biomarkers can depend on observed baseline
    # covariates but contain no latent-age, origin-shift, or frailty information.
    mean_mat <- matrix(0.18 * X1 + 0.10 * X2, nrow = n, ncol = p)
    out <- mean_mat + matrix(
      stats::rnorm(n * p, sd = sc_val(sc, "sigma_bio", SIGMA_BIO_WEAK)),
      nrow = n
    )
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

  beta_U <- sc_val(sc, "risk_beta_U", SURV_BETA_U)
  risk_age_growth <- sc_val(sc, "risk_age_growth", silk_opt("SURV_AGE_GROWTH"))
  eta_risk <- SURV_BETA_X1 * X1 + SURV_BETA_X2 * X2 + beta_U * U_latent

  D_star <- r_attained_age_ph(A_star, eta_risk, growth = risk_age_growth)
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
    risk_age_growth = risk_age_growth,
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
  stop_age <- outer(subjects$A_star, u_grid, "+")
  start_hazard <- attained_age_cumulative_hazard(
    subjects$A_star, subjects$risk_age_growth[1L]
  )
  stop_hazard <- attained_age_cumulative_hazard(
    as.vector(stop_age), subjects$risk_age_growth[1L]
  )
  increment <- matrix(stop_hazard, nrow = nrow(subjects), ncol = length(u_grid)) -
    matrix(start_hazard, nrow = nrow(subjects), ncol = length(u_grid))
  exp(-pmax(increment, 0) * exp(subjects$eta_risk))
}
