# =============================================================================
# registration.R
# Characteristic-kernel SILK registration
# =============================================================================

# The biomarker loss implemented here is exactly
#
#   ||phi(B) - mu_hat(x)||^2_HB,
#
# where k_B is a supported exact characteristic kernel and
# mu_hat(x) = sum_i beta_i(x) phi(B_i).  Its observable expansion is
#
#   1 - 2 beta(x)' k_B(B_train, B) + beta(x)' K_B beta(x).
#
# No engineered biomarker moments and no finite biomarker feature map are used.
# A fixed Fourier map is used only to define the scalar template-input kernel.
# The corresponding vector-valued ridge solution is evaluated in an exactly
# equivalent primal form; the biomarker kernel remains exact and characteristic.

# ----------------------------- Kernel utilities -----------------------------

#' @keywords internal
make_scaler <- function(x) {
  x <- as.matrix(x)
  if (!nrow(x) || !ncol(x) || any(!is.finite(x))) {
    stop("Biomarker matrix must be nonempty and finite.", call. = FALSE)
  }
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  list(center = center, scale = scale)
}

#' @keywords internal
apply_scaler <- function(x, scaler) {
  x <- as.matrix(x)
  if (ncol(x) != length(scaler$center)) {
    stop("Biomarker dimension does not match the fitted scaler.", call. = FALSE)
  }
  x <- sweep(x, 2L, scaler$center, "-")
  x <- sweep(x, 2L, scaler$scale, "/")
  if (any(!is.finite(x))) {
    stop("Scaled biomarkers contain non-finite values.", call. = FALSE)
  }
  x
}

#' @keywords internal
pairwise_squared_distance <- function(x, y = NULL) {
  x <- as.matrix(x)
  if (is.null(y)) y <- x
  y <- as.matrix(y)
  if (ncol(x) != ncol(y)) {
    stop("Kernel inputs have different dimensions.", call. = FALSE)
  }
  distance <- outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y)
  pmax(distance, 0)
}

#' @keywords internal
gaussian_rbf_kernel <- function(x, y = NULL, bandwidth) {
  bandwidth <- as.numeric(bandwidth)[1L]
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    stop("Gaussian biomarker-kernel bandwidth must be strictly positive.", call. = FALSE)
  }
  exp(-pairwise_squared_distance(x, y) / (2 * bandwidth^2))
}

#' @keywords internal
laplace_kernel <- function(x, y = NULL, bandwidth) {
  bandwidth <- as.numeric(bandwidth)[1L]
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    stop("Laplace biomarker-kernel bandwidth must be strictly positive.", call. = FALSE)
  }
  exp(-sqrt(pairwise_squared_distance(x, y)) / bandwidth)
}

#' @keywords internal
matern32_kernel <- function(x, y = NULL, bandwidth) {
  bandwidth <- as.numeric(bandwidth)[1L]
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    stop("Matern-3/2 biomarker-kernel bandwidth must be strictly positive.", call. = FALSE)
  }
  radius <- sqrt(3 * pairwise_squared_distance(x, y)) / bandwidth
  (1 + radius) * exp(-radius)
}

#' @keywords internal
normalize_biomarker_kernel <- function(kernel = NULL) {
  if (is.null(kernel) || !length(kernel) || !nzchar(trimws(as.character(kernel)[1L]))) {
    kernel <- silk_opt("BIOMARKER_KERNEL")
  }
  key <- tolower(gsub("[^a-z0-9]", "", trimws(as.character(kernel)[1L])))
  canonical <- switch(
    key,
    gaussian = "gaussian_rbf",
    gaussianrbf = "gaussian_rbf",
    rbf = "gaussian_rbf",
    laplace = "laplace_l2",
    laplacian = "laplace_l2",
    laplacel2 = "laplace_l2",
    matern12 = "laplace_l2",
    matern = "matern_3_2",
    matern32 = "matern_3_2",
    stop(
      "Unsupported biomarker kernel '", as.character(kernel)[1L],
      "'. Use one of: gaussian, laplace, or matern32. Linear and polynomial ",
      "kernels are deliberately unavailable because they are not generally characteristic.",
      call. = FALSE
    )
  )
  canonical
}

#' @keywords internal
evaluate_biomarker_kernel <- function(x, y = NULL, builder) {
  if (is.null(builder$kernel) || is.null(builder$bandwidth)) {
    stop("Invalid biomarker-kernel builder.", call. = FALSE)
  }
  switch(
    normalize_biomarker_kernel(builder$kernel),
    gaussian_rbf = gaussian_rbf_kernel(x, y, builder$bandwidth),
    laplace_l2 = laplace_kernel(x, y, builder$bandwidth),
    matern_3_2 = matern32_kernel(x, y, builder$bandwidth),
    stop("Unsupported fitted biomarker kernel.", call. = FALSE)
  )
}

#' @keywords internal
biomarker_kernel_label <- function(kernel) {
  switch(
    normalize_biomarker_kernel(kernel),
    gaussian_rbf = "Gaussian RBF",
    laplace_l2 = "Laplace",
    matern_3_2 = "Matern-3/2"
  )
}

#' @keywords internal
median_distance_bandwidth <- function(x, max_points = NULL) {
  x <- as.matrix(x)
  if (is.null(max_points)) max_points <- silk_opt("BIOMARKER_BANDWIDTH_MAX_POINTS")
  n <- nrow(x)
  if (n < 2L) return(1)
  max_points <- suppressWarnings(as.integer(max_points)[1L])
  if (!is.finite(max_points) || max_points < 2L) max_points <- 500L
  m <- min(n, max_points)
  index <- unique(as.integer(round(seq(1, n, length.out = m))))
  distance <- pairwise_squared_distance(x[index, , drop = FALSE])
  values <- distance[upper.tri(distance) & distance > 1e-12 & is.finite(distance)]
  if (!length(values)) return(1)
  bandwidth <- sqrt(stats::median(values))
  if (!is.finite(bandwidth) || bandwidth <= 0) 1 else bandwidth
}

#' @keywords internal
make_biomarker_kernel_builder <- function(visits_df, biomarker_kernel = NULL) {
  columns <- bio_columns(visits_df)
  biomarkers <- as.matrix(visits_df[, columns, drop = FALSE])
  scaler <- make_scaler(biomarkers)
  standardized <- apply_scaler(biomarkers, scaler)
  specification <- silk_opt("BIOMARKER_BANDWIDTH")
  # Multiplier on the selected bandwidth. The confirmatory default is 1, which
  # reproduces the median heuristic exactly. Values other than 1 are used only
  # to run a prespecified bandwidth ladder as a sensitivity analysis; the median
  # heuristic is measured in standardized biomarker-space units, not years, so
  # it is not assumed to be misscaled until the ladder shows it is.
  bandwidth_scale <- suppressWarnings(
    as.numeric(silk_opt("BIOMARKER_BANDWIDTH_SCALE"))[1L]
  )
  if (!is.finite(bandwidth_scale) || bandwidth_scale <= 0) bandwidth_scale <- 1
  numeric_bandwidth <- suppressWarnings(as.numeric(specification)[1L])
  if (is.finite(numeric_bandwidth) && numeric_bandwidth > 0) {
    bandwidth <- numeric_bandwidth
    bandwidth_rule <- "fixed"
  } else if (identical(tolower(trimws(as.character(specification)[1L])), "median")) {
    bandwidth <- median_distance_bandwidth(standardized) * bandwidth_scale
    bandwidth_rule <- if (isTRUE(all.equal(bandwidth_scale, 1))) {
      "median_pairwise_distance"
    } else {
      paste0("median_pairwise_distance_x", format(bandwidth_scale))
    }
  } else {
    stop(
      "BIOMARKER_BANDWIDTH must be a positive number or 'median'.",
      call. = FALSE
    )
  }
  list(
    kernel = normalize_biomarker_kernel(biomarker_kernel),
    characteristic = TRUE,
    normalized_diagonal = TRUE,
    biomarker_cols = columns,
    scaler = scaler,
    bandwidth = bandwidth,
    bandwidth_rule = bandwidth_rule
  )
}

#' @keywords internal
apply_biomarker_builder <- function(visits_df, builder) {
  missing_columns <- setdiff(builder$biomarker_cols, names(visits_df))
  if (length(missing_columns)) {
    stop(
      "Missing biomarker columns: ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  apply_scaler(
    as.matrix(visits_df[, builder$biomarker_cols, drop = FALSE]),
    builder$scaler
  )
}

#' @keywords internal
longitudinal_kernel_signal <- function(visits_df, standardized, builder) {
  similarities <- numeric(0)
  negative_lag_distance <- numeric(0)
  for (subject_id in unique(visits_df$id)) {
    index <- which(visits_df$id == subject_id)
    if (length(index) < 2L) next
    kernel <- evaluate_biomarker_kernel(
      standardized[index, , drop = FALSE], builder = builder
    )
    lag_distance <- abs(outer(visits_df$lag[index], visits_df$lag[index], "-"))
    upper <- upper.tri(kernel)
    similarities <- c(similarities, kernel[upper])
    negative_lag_distance <- c(negative_lag_distance, -lag_distance[upper])
  }
  # A zero-variance lag-gap design (for example a balanced two-visit schedule,
  # in which every subject contributes exactly one pair at the same lag gap)
  # makes this correlation mathematically undefined. That is a property of the
  # visit design, not evidence that biomarkers carry no longitudinal signal, so
  # the statistic returns NA_real_ ("not evaluable") rather than 0. Callers must
  # treat NA as "gate not evaluable, do not disable"; see the guards in
  # biomarker_clock_initial_shift() and template_clock_initial_shift().
  if (length(similarities) < 3L || stats::sd(similarities) < 1e-12 ||
      stats::sd(negative_lag_distance) < 1e-12) return(NA_real_)
  value <- suppressWarnings(stats::cor(similarities, negative_lag_distance))
  if (is.finite(value)) value else NA_real_
}

#' Should the longitudinal-signal gate disable the biomarker clock?
#'
#' Returns TRUE only when the gate statistic is evaluable and falls below
#' CLOCK_SIGNAL_MIN. A non-evaluable statistic (NA) leaves the clock enabled:
#' an undefined diagnostic is not evidence of absent signal.
#'
#' @param longitudinal_signal Numeric scalar or NA from
#'   \code{longitudinal_kernel_signal}.
#' @return Logical scalar.
#' @keywords internal
clock_gate_disables <- function(longitudinal_signal) {
  minimum <- suppressWarnings(as.numeric(silk_opt("CLOCK_SIGNAL_MIN")))
  if (!length(minimum) || !is.finite(minimum[1L])) return(FALSE)
  # Length-zero inputs (NULL, numeric(0), a dropped list element) are as
  # non-evaluable as NA and must not be fed to `||`, which errors on
  # zero-length operands. Reduce to a single finite value or bail out.
  if (!length(longitudinal_signal)) return(FALSE)
  value <- suppressWarnings(as.numeric(longitudinal_signal)[1L])
  if (!is.finite(value)) return(FALSE)
  value < minimum[1L]
}

#' @keywords internal
with_preserved_seed <- function(seed, expr) {
  if (is.null(seed) || !length(seed) || !is.finite(seed[1L])) return(force(expr))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed[1L]))
  force(expr)
}

#' @keywords internal
template_covariate_bandwidth <- function(name) {
  value <- switch(
    name,
    X1 = silk_opt("H_X_BANDWIDTH"),
    X2 = silk_opt("H_X2_BANDWIDTH"),
    X3 = silk_opt("H_X3_BANDWIDTH"),
    X4 = silk_opt("H_X4_BANDWIDTH"),
    lag = silk_opt("H_LAG_BANDWIDTH"),
    1
  )
  value <- as.numeric(value)[1L]
  if (!is.finite(value) || value <= 0) {
    stop("Template-input bandwidth for ", name, " must be positive.", call. = FALSE)
  }
  value
}

#' @keywords internal
make_input_kernel_builder <- function(visits_df) {
  dimension <- suppressWarnings(as.integer(silk_opt("TEMPLATE_INPUT_FEATURES"))[1L])
  if (!is.finite(dimension) || dimension < 2L || dimension %% 2L != 0L) {
    stop("TEMPLATE_INPUT_FEATURES must be an even integer of at least two.", call. = FALSE)
  }
  requested <- as.character(silk_opt("TEMPLATE_INPUT_COVARIATES"))
  covariates <- intersect(requested, names(visits_df))
  input_names <- c("candidate_age", covariates)
  bandwidths <- c(
    candidate_age = as.numeric(silk_opt("H_A_BANDWIDTH"))[1L],
    stats::setNames(vapply(covariates, template_covariate_bandwidth, numeric(1)), covariates)
  )
  if (any(!is.finite(bandwidths)) || any(bandwidths <= 0)) {
    stop("All template-input bandwidths must be positive and finite.", call. = FALSE)
  }
  n_frequency <- dimension %/% 2L
  frequency <- with_preserved_seed(
    silk_opt("TEMPLATE_FEATURE_SEED"),
    matrix(
      stats::rnorm(length(input_names) * n_frequency),
      nrow = length(input_names), ncol = n_frequency
    )
  )
  list(
    kernel = "fixed_fourier_input_kernel",
    dimension = dimension,
    n_frequency = n_frequency,
    frequency = frequency,
    covariates = covariates,
    input_names = input_names,
    bandwidths = bandwidths
  )
}

#' @keywords internal
template_input_matrix <- function(candidate_age, visits_df, builder) {
  candidate_age <- as.numeric(candidate_age)
  if (length(candidate_age) != nrow(visits_df)) {
    stop("candidate_age must have one value per visit row.", call. = FALSE)
  }
  if (length(builder$covariates)) {
    missing_columns <- setdiff(builder$covariates, names(visits_df))
    if (length(missing_columns)) {
      stop(
        "Missing template-input columns: ", paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    input <- cbind(
      candidate_age = candidate_age,
      as.matrix(visits_df[, builder$covariates, drop = FALSE])
    )
  } else {
    input <- matrix(candidate_age, ncol = 1L, dimnames = list(NULL, "candidate_age"))
  }
  storage.mode(input) <- "double"
  input <- sweep(input, 2L, builder$bandwidths, "/")
  if (any(!is.finite(input))) {
    stop("Template inputs contain non-finite values.", call. = FALSE)
  }
  input
}

#' @keywords internal
apply_input_kernel_builder <- function(candidate_age, visits_df, builder) {
  input <- template_input_matrix(candidate_age, visits_df, builder)
  angle <- input %*% builder$frequency
  sqrt(1 / builder$n_frequency) * cbind(cos(angle), sin(angle))
}

#' @keywords internal
safe_chol <- function(matrix, base_jitter = NULL) {
  if (is.null(base_jitter)) base_jitter <- silk_opt("TEMPLATE_NUMERICAL_JITTER")
  matrix <- (matrix + t(matrix)) / 2
  jitter <- max(0, as.numeric(base_jitter)[1L])
  for (attempt in 0:7) {
    addition <- if (attempt == 0L) 0 else max(jitter, 1e-12) * 10^(attempt - 1L)
    factor <- tryCatch(
      chol(matrix + addition * diag(nrow(matrix))),
      error = function(e) NULL
    )
    if (!is.null(factor)) return(list(R = factor, jitter = addition))
  }
  stop("Template ridge system is not numerically positive definite.", call. = FALSE)
}

#' @keywords internal
chol_solve <- function(chol_factor, right_hand_side) {
  backsolve(chol_factor, forwardsolve(t(chol_factor), right_hand_side))
}

# ----------------------------- Parallel helpers -----------------------------

#' @keywords internal
silk_parallel_lapply <- function(X, FUN, ...) {
  n_cores <- min(silk_resolve_cores(), length(X))
  nested <- nzchar(Sys.getenv("SILK_PARALLEL_ACTIVE"))
  if (n_cores <= 1L || nested || !requireNamespace("parallel", quietly = TRUE)) {
    return(lapply(X, FUN, ...))
  }

  Sys.setenv(SILK_PARALLEL_ACTIVE = "1")
  on.exit(Sys.unsetenv("SILK_PARALLEL_ACTIVE"), add = TRUE)

  if (.Platform$OS.type == "unix") {
    result <- parallel::mclapply(
      X, FUN, ..., mc.cores = n_cores, mc.preschedule = FALSE, mc.set.seed = TRUE
    )
    failed <- vapply(result, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      for (index in which(failed)) result[[index]] <- FUN(X[[index]], ...)
    }
    return(result)
  }

  if (!requireNamespace("SILK", quietly = TRUE)) return(lapply(X, FUN, ...))
  cluster <- tryCatch(parallel::makePSOCKcluster(n_cores), error = function(e) NULL)
  if (is.null(cluster)) return(lapply(X, FUN, ...))
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  current_options <- options()[grep("^silk\\.", names(options()), value = TRUE)]
  initialized <- tryCatch({
    parallel::clusterCall(cluster, function(option_values) {
      suppressWarnings(suppressMessages(library(SILK)))
      options(option_values)
      Sys.setenv(SILK_PARALLEL_ACTIVE = "1")
      invisible(NULL)
    }, current_options)
    TRUE
  }, error = function(e) FALSE)
  if (!initialized) return(lapply(X, FUN, ...))
  result <- tryCatch(parallel::parLapply(cluster, X, FUN, ...), error = function(e) NULL)
  if (is.null(result)) lapply(X, FUN, ...) else result
}

# ----------------------------- Anchor helpers -------------------------------

#' @keywords internal
anchor_basis <- function(subjects, ids, mode = NULL) {
  if (is.null(mode)) mode <- silk_opt("ANCHOR_MODE")
  matched <- subjects[match(ids, subjects$id), , drop = FALSE]
  if (anyNA(matched$id)) stop("Registration ids are missing from subjects.", call. = FALSE)
  if (identical(mode, "intercept")) {
    basis <- matrix(1, nrow = length(ids), ncol = 1L)
    colnames(basis) <- "intercept"
    return(basis)
  }
  if (!"X1" %in% names(matched)) {
    stop("ANCHOR_MODE='rx' requires subject-level X1.", call. = FALSE)
  }
  x1 <- as.numeric(scale(matched$X1))
  x1[!is.finite(x1)] <- 0
  cbind(intercept = 1, X1 = x1)
}

#' @keywords internal
project_to_anchor <- function(shift, anchor, lower, upper) {
  anchor <- as.matrix(anchor)
  if (ncol(anchor) == 0L) return(pmin(pmax(shift, lower), upper))
  gram <- crossprod(anchor)
  right_hand_side <- colSums(anchor * shift)
  adjustment <- tryCatch(
    anchor %*% solve(gram + 1e-8 * diag(ncol(anchor)), right_hand_side),
    error = function(e) 0
  )
  pmin(pmax(as.numeric(shift - adjustment), lower), upper)
}

#' @keywords internal
shift_ridge_matrix <- function(grid, n) {
  penalty <- suppressWarnings(as.numeric(silk_opt("SHIFT_RIDGE"))[1L])
  if (!is.finite(penalty) || penalty <= 0) return(NULL)
  matrix(penalty * grid^2, nrow = n, ncol = length(grid), byrow = TRUE)
}

#' @keywords internal
constrained_grid_update <- function(loss, grid, anchor, center_shift = NULL) {
  n <- nrow(loss)
  anchor <- as.matrix(anchor)
  if (nrow(anchor) != n) stop("anchor has the wrong number of rows.", call. = FALSE)
  lower <- min(grid)
  upper <- max(grid)
  step <- as.numeric(silk_opt("SHIFT_GRID_STEP"))
  ridge <- shift_ridge_matrix(grid, n)
  local_radius <- as.numeric(silk_opt("ALT_LOCAL_RADIUS"))[1L]
  if (!is.null(center_shift) && is.finite(local_radius) && local_radius > 0) {
    outside <- abs(outer(as.numeric(center_shift), grid, "-")) > local_radius + 1e-12
    loss[outside] <- Inf
  }

  choose_for_lambda <- function(lambda) {
    score <- as.vector(anchor %*% lambda)
    adjusted <- loss + outer(score, grid, "*")
    if (!is.null(ridge)) adjusted <- adjusted + ridge
    grid[max.col(-adjusted, ties.method = "first")]
  }

  if (ncol(anchor) == 1L) {
    moment <- function(lambda) sum(anchor[, 1L] * choose_for_lambda(lambda))
    lower_lambda <- -1
    upper_lambda <- 1
    lower_moment <- moment(lower_lambda)
    upper_moment <- moment(upper_lambda)
    expansion <- 0L
    while (is.finite(lower_moment) && is.finite(upper_moment) &&
           lower_moment * upper_moment > 0 && expansion < 40L) {
      lower_lambda <- lower_lambda * 2
      upper_lambda <- upper_lambda * 2
      lower_moment <- moment(lower_lambda)
      upper_moment <- moment(upper_lambda)
      expansion <- expansion + 1L
    }
    if (is.finite(lower_moment) && is.finite(upper_moment) &&
        lower_moment * upper_moment <= 0) {
      for (iteration in seq_len(60L)) {
        midpoint <- (lower_lambda + upper_lambda) / 2
        midpoint_moment <- moment(midpoint)
        if (abs(midpoint_moment) <= max(1e-8, 0.01 * n * step)) break
        if (lower_moment * midpoint_moment <= 0) {
          upper_lambda <- midpoint
          upper_moment <- midpoint_moment
        } else {
          lower_lambda <- midpoint
          lower_moment <- midpoint_moment
        }
      }
      shift <- choose_for_lambda((lower_lambda + upper_lambda) / 2)
    } else {
      selection_loss <- if (is.null(ridge)) loss else loss + ridge
      shift <- grid[apply(selection_loss, 1L, which.min)]
    }
    return(project_to_anchor(shift, anchor, lower, upper))
  }

  objective <- function(lambda) {
    shift <- choose_for_lambda(lambda)
    moment <- colSums(anchor * shift)
    sum(moment^2) + 1e-8 * sum(lambda^2)
  }
  starts <- list(rep(0, ncol(anchor)))
  for (index in seq_len(4L)) {
    starts[[length(starts) + 1L]] <- stats::rnorm(ncol(anchor), sd = 0.5 * index)
  }
  best <- NULL
  best_value <- Inf
  for (start in starts) {
    fit <- tryCatch(
      stats::optim(start, objective, method = "Nelder-Mead", control = list(maxit = 200)),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$value) && fit$value < best_value) {
      best <- fit$par
      best_value <- fit$value
    }
  }
  if (is.null(best)) best <- rep(0, ncol(anchor))
  project_to_anchor(choose_for_lambda(best), anchor, lower, upper)
}

# ----------------------------- RKHS template --------------------------------

#' @keywords internal
prepare_registration_data <- function(visits_df, biomarker_kernel = NULL) {
  required <- c("id", "visit", "A_obs_il")
  missing_columns <- setdiff(required, names(visits_df))
  if (length(missing_columns)) {
    stop("Missing visit columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  biomarker_builder <- make_biomarker_kernel_builder(visits_df, biomarker_kernel)
  standardized <- apply_biomarker_builder(visits_df, biomarker_builder)
  input_builder <- make_input_kernel_builder(visits_df)
  positions <- sort(unique(visits_df$visit))
  by_position <- stats::setNames(vector("list", length(positions)), as.character(positions))
  for (position in positions) {
    index <- which(visits_df$visit == position)
    biomarkers <- standardized[index, , drop = FALSE]
    by_position[[as.character(position)]] <- list(
      row_index = index,
      B_std = biomarkers,
      K_B = evaluate_biomarker_kernel(biomarkers, builder = biomarker_builder)
    )
  }
  list(
    biomarker_builder = biomarker_builder,
    input_builder = input_builder,
    B_std = standardized,
    longitudinal_signal = longitudinal_kernel_signal(
      visits_df, standardized, biomarker_builder
    ),
    positions = positions,
    by_position = by_position
  )
}

#' @keywords internal
build_position_template <- function(visits_df, subjects, shift, ids,
                                    registration_data = NULL) {
  if (is.null(registration_data)) registration_data <- prepare_registration_data(visits_df)
  id_index <- match(visits_df$id, ids)
  if (anyNA(id_index)) stop("Training visits contain an id absent from ids.", call. = FALSE)
  calibrated_age <- visits_df$A_obs_il - shift[id_index]
  number_positions <- length(registration_data$positions)
  lambda <- as.numeric(silk_opt("TEMPLATE_RIDGE_LAMBDA"))[1L]
  if (!is.finite(lambda) || lambda <= 0) {
    stop("TEMPLATE_RIDGE_LAMBDA must be positive and finite.", call. = FALSE)
  }
  templates <- list()
  clock_sse <- 0
  clock_sst <- 0

  for (position in registration_data$positions) {
    key <- as.character(position)
    cached <- registration_data$by_position[[key]]
    index <- cached$row_index
    position_visits <- visits_df[index, , drop = FALSE]
    psi <- apply_input_kernel_builder(
      calibrated_age[index], position_visits, registration_data$input_builder
    )
    rho <- length(ids) * number_positions * lambda
    ridge_system <- crossprod(psi) + rho * diag(ncol(psi))
    factor <- safe_chol(ridge_system)

    # A = Psi(Psi'Psi + rho I)^(-1). For K_X = Psi Psi',
    # (K_X + rho I)^(-1) k_X(x) = A psi(x) exactly.
    coefficient <- chol_solve(factor$R, t(psi))
    A <- t(coefficient)
    biomarker_quadratic <- crossprod(A, cached$K_B %*% A)
    biomarker_quadratic <- (biomarker_quadratic + t(biomarker_quadratic)) / 2

    templates[[key]] <- list(
      visit = position,
      B_train = cached$B_std,
      landmark_observed = subjects$A_obs[match(position_visits$id, subjects$id)],
      K_B = cached$K_B,
      A = A,
      M_B = biomarker_quadratic,
      penalty_norm = sum(diag(biomarker_quadratic)),
      rho = rho,
      numerical_jitter = factor$jitter
    )
    clock_response <- templates[[key]]$landmark_observed
    clock_prediction <- kernel_neighbour_prediction(
      cached$K_B, clock_response, exclude_self = TRUE
    )
    clock_sse <- clock_sse + sum((clock_response - clock_prediction)^2)
    clock_sst <- clock_sst + sum((clock_response - mean(clock_response))^2)
  }

  list(
    templates = templates,
    ids = ids,
    positions = registration_data$positions,
    biomarker_builder = registration_data$biomarker_builder,
    input_builder = registration_data$input_builder,
    biomarker_kernel = registration_data$biomarker_builder$kernel,
    characteristic = TRUE,
    clock_signal_r2 = if (clock_sst > 0) 1 - clock_sse / clock_sst else NA_real_,
    longitudinal_signal = registration_data$longitudinal_signal,
    template_ridge = lambda,
    number_positions = number_positions
  )
}

#' @keywords internal
position_cross_basis <- function(template_position, query_biomarkers, biomarker_builder) {
  same_rows <- nrow(query_biomarkers) == nrow(template_position$B_train) &&
    ncol(query_biomarkers) == ncol(template_position$B_train) &&
    isTRUE(all.equal(
      query_biomarkers, template_position$B_train,
      tolerance = 0, check.attributes = FALSE
    ))
  cross_kernel <- if (same_rows) {
    template_position$K_B
  } else {
    evaluate_biomarker_kernel(
      template_position$B_train, query_biomarkers, builder = biomarker_builder
    )
  }
  crossprod(template_position$A, cross_kernel)
}

#' @keywords internal
rkhs_loss_from_summaries <- function(psi, cross_summary, biomarker_quadratic) {
  cross_term <- rowSums(psi * cross_summary)
  template_norm <- rowSums(psi * (psi %*% biomarker_quadratic))
  # k_B(b,b)=1 for every supported normalized kernel. Negative values below numerical zero
  # arise only from floating-point error in a positive-semidefinite expression.
  pmax(1 - 2 * cross_term + template_norm, 0)
}

#' @keywords internal
position_loss_grid <- function(template_position, visits, query_biomarkers, grid,
                               input_builder, biomarker_builder) {
  number_rows <- nrow(visits)
  number_grid <- length(grid)
  output <- matrix(NA_real_, nrow = number_rows, ncol = number_grid)
  cross_basis <- position_cross_basis(template_position, query_biomarkers, biomarker_builder)
  query_chunk <- suppressWarnings(as.integer(silk_opt("TEMPLATE_QUERY_CHUNK"))[1L])
  if (!is.finite(query_chunk) || query_chunk < 1L) query_chunk <- 5000L
  grids_per_chunk <- max(1L, floor(query_chunk / max(1L, number_rows)))

  for (start in seq(1L, number_grid, by = grids_per_chunk)) {
    end <- min(number_grid, start + grids_per_chunk - 1L)
    grid_block <- grid[start:end]
    repeat_index <- rep(seq_len(number_rows), times = length(grid_block))
    query_visits <- visits[repeat_index, , drop = FALSE]
    query_age <- as.vector(outer(visits$A_obs_il, grid_block, "-"))
    psi <- apply_input_kernel_builder(query_age, query_visits, input_builder)
    cross_summary <- t(cross_basis[, repeat_index, drop = FALSE])
    loss <- rkhs_loss_from_summaries(psi, cross_summary, template_position$M_B)
    output[, start:end] <- matrix(loss, nrow = number_rows, ncol = length(grid_block))
  }
  output
}

#' @keywords internal
position_loss_vector <- function(template_position, visits, query_biomarkers, shift,
                                 input_builder, biomarker_builder) {
  cross_basis <- position_cross_basis(template_position, query_biomarkers, biomarker_builder)
  psi <- apply_input_kernel_builder(visits$A_obs_il - shift, visits, input_builder)
  rkhs_loss_from_summaries(psi, t(cross_basis), template_position$M_B)
}

#' @keywords internal
loss_grid_from_template <- function(template, visits_df, grid) {
  ids <- sort(unique(visits_df$id))
  number_ids <- length(ids)
  standardized <- apply_biomarker_builder(visits_df, template$biomarker_builder)
  total_loss <- matrix(0, nrow = number_ids, ncol = length(grid))
  counts <- matrix(0, nrow = number_ids, ncol = length(grid))

  for (position in sort(unique(visits_df$visit))) {
    index <- which(visits_df$visit == position)
    template_position <- template$templates[[as.character(position)]]
    if (!length(index) || is.null(template_position)) next
    position_visits <- visits_df[index, , drop = FALSE]
    row_id <- match(position_visits$id, ids)
    loss <- position_loss_grid(
      template_position, position_visits, standardized[index, , drop = FALSE], grid,
      template$input_builder, template$biomarker_builder
    )
    aggregated <- rowsum(loss, row_id, reorder = FALSE)
    rows <- as.integer(rownames(aggregated))
    row_counts <- tabulate(row_id, nbins = number_ids)[rows]
    total_loss[rows, ] <- total_loss[rows, , drop = FALSE] + aggregated
    counts[rows, ] <- counts[rows, , drop = FALSE] + row_counts
  }

  if (any(rowSums(counts) == 0)) {
    stop("At least one subject has no visit position represented in the template.", call. = FALSE)
  }
  total_loss / pmax(counts, 1)
}

#' @keywords internal
registration_total_objective <- function(template, visits_df, shift, grid_dummy = NULL) {
  ids <- sort(unique(visits_df$id))
  if (length(shift) != length(ids)) stop("shift has the wrong length.", call. = FALSE)
  standardized <- apply_biomarker_builder(visits_df, template$biomarker_builder)
  loss <- numeric(length(ids))
  count <- numeric(length(ids))

  for (position in sort(unique(visits_df$visit))) {
    index <- which(visits_df$visit == position)
    template_position <- template$templates[[as.character(position)]]
    if (!length(index) || is.null(template_position)) next
    position_visits <- visits_df[index, , drop = FALSE]
    row_id <- match(position_visits$id, ids)
    point_loss <- position_loss_vector(
      template_position, position_visits, standardized[index, , drop = FALSE],
      shift[row_id], template$input_builder, template$biomarker_builder
    )
    aggregated <- rowsum(point_loss, row_id, reorder = FALSE)
    rows <- as.integer(rownames(aggregated))
    loss[rows] <- loss[rows] + as.vector(aggregated)
    count <- count + tabulate(row_id, nbins = length(ids))
  }

  data_term <- mean(loss / pmax(count, 1))
  penalty <- template$template_ridge * sum(vapply(
    template$templates, function(value) value$penalty_norm, numeric(1)
  ))
  data_term + penalty
}

# ----------------------------- Alternating fit ------------------------------

#' @keywords internal
fit_registration_train <- function(visits_df, subjects, init_e = NULL, grid,
                                   max_iter = NULL, tol = NULL,
                                   registration_data = NULL) {
  if (is.null(max_iter)) max_iter <- silk_opt("N_ALT_ITER_MAX")
  if (is.null(tol)) tol <- silk_opt("ALT_TOL")
  ids <- sort(unique(visits_df$id))
  number_ids <- length(ids)
  if (is.null(init_e)) init_e <- rep(0, number_ids)
  if (length(init_e) != number_ids) stop("init_e has the wrong length.", call. = FALSE)
  if (is.null(registration_data)) registration_data <- prepare_registration_data(visits_df)
  anchor <- anchor_basis(subjects, ids)
  shift <- project_to_anchor(as.numeric(init_e), anchor, min(grid), max(grid))
  objective_trace <- numeric(0)
  step <- as.numeric(silk_opt("SHIFT_GRID_STEP"))

  for (iteration in seq_len(max_iter)) {
    template <- build_position_template(
      visits_df, subjects, shift, ids, registration_data = registration_data
    )
    loss <- loss_grid_from_template(template, visits_df, grid)
    new_shift <- constrained_grid_update(loss, grid, anchor, center_shift = shift)
    objective <- registration_total_objective(template, visits_df, new_shift)
    objective_trace <- c(objective_trace, objective)
    max_change <- max(abs(new_shift - shift))
    shift <- new_shift
    if (length(objective_trace) >= 2L &&
        abs(diff(utils::tail(objective_trace, 2L))) < tol && max_change < step) break
  }

  template <- build_position_template(
    visits_df, subjects, shift, ids, registration_data = registration_data
  )
  list(
    e_train = shift,
    template = template,
    ids = ids,
    obj_trace = objective_trace,
    grid = grid,
    implementation = paste0(
      "exact characteristic ",
      biomarker_kernel_label(template$biomarker_builder$kernel),
      " biomarker-kernel RKHS loss"
    )
  )
}

#' @keywords internal
kernel_neighbour_prediction <- function(kernel, response, exclude_self = FALSE) {
  kernel <- as.matrix(kernel)
  response <- as.numeric(response)
  if (ncol(kernel) != length(response)) {
    stop("Kernel columns and clock responses have different lengths.", call. = FALSE)
  }
  if (isTRUE(exclude_self) && nrow(kernel) == ncol(kernel)) diag(kernel) <- 0
  available <- ncol(kernel) - as.integer(isTRUE(exclude_self) && nrow(kernel) == ncol(kernel))
  number_neighbours <- min(available, max(8L, ceiling(sqrt(ncol(kernel)))))
  if (number_neighbours < 1L) return(rep(mean(response), nrow(kernel)))
  vapply(seq_len(nrow(kernel)), function(row) {
    order_index <- order(kernel[row, ], decreasing = TRUE)[seq_len(number_neighbours)]
    weight <- pmax(kernel[row, order_index], 0)
    if (!any(is.finite(weight)) || sum(weight, na.rm = TRUE) <= 1e-12) {
      return(mean(response))
    }
    sum(weight * response[order_index]) / sum(weight)
  }, numeric(1))
}

#' @keywords internal
biomarker_clock_initial_shift <- function(visits_df, subjects, ids,
                                          registration_data, grid) {
  subject_landmark <- subjects$A_obs[match(ids, subjects$id)]
  stage_sum <- numeric(length(ids))
  stage_count <- integer(length(ids))
  clock_sse <- 0
  clock_sst <- 0
  for (position in registration_data$positions) {
    cached <- registration_data$by_position[[as.character(position)]]
    index <- cached$row_index
    if (length(index) < 3L) next
    position_visits <- visits_df[index, , drop = FALSE]
    row_subject <- match(position_visits$id, ids)
    response <- subject_landmark[row_subject]
    predicted_landmark <- kernel_neighbour_prediction(
      cached$K_B, response, exclude_self = TRUE
    )
    clock_sse <- clock_sse + sum((response - predicted_landmark)^2)
    clock_sst <- clock_sst + sum((response - mean(response))^2)
    stage_sum[row_subject] <- stage_sum[row_subject] + predicted_landmark
    stage_count[row_subject] <- stage_count[row_subject] + 1L
  }
  fallback <- mean(subject_landmark)
  stage <- ifelse(stage_count > 0L, stage_sum / pmax(stage_count, 1L), fallback)
  signal_r2 <- if (clock_sst > 0) 1 - clock_sse / clock_sst else NA_real_
  longitudinal_signal <- registration_data$longitudinal_signal
  if (clock_gate_disables(longitudinal_signal)) {
    stage <- subject_landmark
  }
  anchor <- anchor_basis(subjects, ids)
  shift <- project_to_anchor(subject_landmark - stage, anchor, min(grid), max(grid))
  attr(shift, "clock_signal_r2") <- signal_r2
  shift
}

#' @keywords internal
fit_registration_multistart <- function(visits_df, subjects, grid, seed = NULL,
                                        biomarker_kernel = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ids <- sort(unique(visits_df$id))
  number_ids <- length(ids)
  anchor <- anchor_basis(subjects, ids)
  registration_data <- prepare_registration_data(visits_df, biomarker_kernel)
  number_starts <- suppressWarnings(as.integer(silk_opt("N_STARTS"))[1L])
  if (!is.finite(number_starts) || number_starts < 1L) number_starts <- 1L
  random_start_sd <- as.numeric(silk_opt("RANDOM_START_SD"))[1L]

  starts <- list(
    biology_clock = biomarker_clock_initial_shift(
      visits_df, subjects, ids, registration_data, grid
    )
  )
  if (number_starts >= 2L) starts$zero <- rep(0, number_ids)
  if (number_starts >= 3L) {
    for (index in 3:number_starts) {
      starts[[paste0("random", index - 2L)]] <- project_to_anchor(
        stats::rnorm(number_ids, sd = (index - 2L) * random_start_sd),
        anchor, min(grid), max(grid)
      )
    }
  }

  best <- NULL
  best_objective <- Inf
  best_start <- NA_character_
  # The manuscript's diagnostics section promises "the spread of attained
  # objectives across starts". Retain every attained objective, not only the
  # winner, so multistart stability can be reported and audited.
  start_objectives <- stats::setNames(
    rep(NA_real_, length(starts)), names(starts)
  )
  for (name in names(starts)) {
    fit <- tryCatch(
      fit_registration_train(
        visits_df, subjects, starts[[name]], grid,
        registration_data = registration_data
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    objective <- registration_total_objective(fit$template, visits_df, fit$e_train)
    start_objectives[[name]] <- objective
    if (is.finite(objective) && objective < best_objective) {
      best <- fit
      best_objective <- objective
      best_start <- name
    }
  }
  if (is.null(best)) stop("Every SILK registration start failed.", call. = FALSE)

  finite_objectives <- start_objectives[is.finite(start_objectives)]
  sorted <- sort(finite_objectives)
  best$best_obj <- best_objective
  best$best_start <- best_start
  best$start_objectives <- start_objectives
  best$multistart <- list(
    n_starts = length(starts),
    n_successful_starts = length(finite_objectives),
    best_start = best_start,
    best_objective = best_objective,
    objective_spread = if (length(finite_objectives) >= 2L) {
      max(finite_objectives) - min(finite_objectives)
    } else NA_real_,
    second_best_gap = if (length(sorted) >= 2L) {
      unname(sorted[2L] - sorted[1L])
    } else NA_real_,
    start_names = paste(names(start_objectives), collapse = "|"),
    start_objective_values = paste(
      formatC(start_objectives, format = "g", digits = 10), collapse = "|"
    )
  )
  best
}

#' @keywords internal
refine_profile_min <- function(loss_row, grid) {
  step <- as.numeric(silk_opt("SHIFT_GRID_STEP"))
  temperature <- as.numeric(silk_opt("PROFILE_TEMPERATURE"))
  radius <- as.numeric(silk_opt("PROFILE_LOCAL_RADIUS"))
  minimum <- which.min(loss_row)
  estimate <- grid[minimum]
  window <- which(abs(grid - estimate) <= max(2 * step, radius))
  if (length(window) >= 3L && length(unique(grid[window])) >= 3L) {
    fit <- tryCatch(
      stats::lm(loss_row[window] ~ grid[window] + I(grid[window]^2)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      coefficient <- stats::coef(fit)
      if (length(coefficient) == 3L && is.finite(coefficient[2L]) &&
          is.finite(coefficient[3L]) && coefficient[3L] > 0) {
        vertex <- -coefficient[2L] / (2 * coefficient[3L])
        if (is.finite(vertex) && vertex >= min(grid[window]) && vertex <= max(grid[window])) {
          estimate <- vertex
        }
      }
    }
  }

  if (is.finite(temperature) && temperature > 0) {
    local <- which(abs(grid - grid[minimum]) <= radius)
    relative <- loss_row[local] - min(loss_row[local], na.rm = TRUE)
    if (is.finite(max(relative, na.rm = TRUE)) &&
        max(relative, na.rm = TRUE) > 2 * temperature) {
      weight <- exp(-relative / temperature)
      if (sum(weight) > 0 && all(is.finite(weight))) {
        estimate <- sum(grid[local] * weight) / sum(weight)
      }
    }
  }
  pmin(pmax(estimate, min(grid)), max(grid))
}

#' @keywords internal
template_clock_initial_shift <- function(template, visits_df) {
  ids <- sort(unique(visits_df$id))
  if (clock_gate_disables(template$longitudinal_signal)) {
    return(stats::setNames(rep(0, length(ids)), ids))
  }
  standardized <- apply_biomarker_builder(visits_df, template$biomarker_builder)
  stage_sum <- numeric(length(ids))
  stage_count <- integer(length(ids))
  observed_landmark_sum <- numeric(length(ids))

  for (position in sort(unique(visits_df$visit))) {
    index <- which(visits_df$visit == position)
    template_position <- template$templates[[as.character(position)]]
    if (!length(index) || is.null(template_position)) next
    position_visits <- visits_df[index, , drop = FALSE]
    row_id <- match(position_visits$id, ids)
    cross_kernel <- evaluate_biomarker_kernel(
      standardized[index, , drop = FALSE], template_position$B_train,
      builder = template$biomarker_builder
    )
    predicted_landmark <- kernel_neighbour_prediction(
      cross_kernel, template_position$landmark_observed, exclude_self = FALSE
    )
    stage_sum[row_id] <- stage_sum[row_id] + predicted_landmark
    observed_landmark_sum[row_id] <- observed_landmark_sum[row_id] +
      position_visits$A_obs_il + position_visits$lag
    stage_count[row_id] <- stage_count[row_id] + 1L
  }
  observed_landmark <- observed_landmark_sum / pmax(stage_count, 1L)
  stage <- stage_sum / pmax(stage_count, 1L)
  shift <- observed_landmark - stage
  shift[stage_count == 0L | !is.finite(shift)] <- 0
  stats::setNames(shift, ids)
}

#' @keywords internal
predict_registration_shift <- function(template, visits_df, grid) {
  ids <- sort(unique(visits_df$id))
  loss <- loss_grid_from_template(template, visits_df, grid)
  ridge <- shift_ridge_matrix(grid, nrow(loss))
  selection_loss <- if (is.null(ridge)) loss else loss + ridge
  clock_shift <- template_clock_initial_shift(template, visits_df)
  clock_center <- pmin(pmax(
    as.numeric(clock_shift[match(ids, names(clock_shift))]), min(grid)
  ), max(grid))
  search_radius <- as.numeric(silk_opt("PROFILE_SEARCH_RADIUS"))[1L]
  if (is.finite(search_radius) && search_radius > 0) {
    outside <- abs(outer(clock_center, grid, "-")) >
      search_radius + 1e-12
    selection_loss[outside] <- Inf
  }
  best <- max.col(-selection_loss, ties.method = "first")
  refined <- vapply(
    seq_along(ids), function(index) refine_profile_min(selection_loss[index, ], grid),
    numeric(1)
  )
  data.frame(
    id = ids,
    e_hat = refined,
    e_clock = clock_center,
    e_hat_grid = grid[best],
    min_loss = loss[cbind(seq_along(ids), best)],
    at_boundary = best %in% c(1L, length(grid)),
    gap_q1 = vapply(seq_along(ids), function(index) {
      far <- abs(grid - grid[best[index]]) >= 1
      if (!any(far)) return(NA_real_)
      min(loss[index, far]) - loss[index, best[index]]
    }, numeric(1)),
    gap_q2 = vapply(seq_along(ids), function(index) {
      far <- abs(grid - grid[best[index]]) >= 2
      if (!any(far)) return(NA_real_)
      min(loss[index, far]) - loss[index, best[index]]
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
crossfit_registration <- function(train_visits, train_subjects, grid, seed = NULL,
                                  biomarker_kernel = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ids <- sort(unique(train_visits$id))
  number_ids <- length(ids)
  if (number_ids < 2L) stop("Cross-fitted registration needs at least two subjects.", call. = FALSE)
  requested_folds <- suppressWarnings(as.integer(silk_opt("N_FOLDS"))[1L])
  if (!is.finite(requested_folds) || requested_folds < 2L) requested_folds <- 2L
  number_folds <- min(requested_folds, number_ids)
  biomarker_kernel <- normalize_biomarker_kernel(biomarker_kernel)
  fold <- sample(rep(seq_len(number_folds), length.out = number_ids))
  output <- data.frame(
    id = ids, e_hat = NA_real_, e_clock = NA_real_, S_hat = NA_real_, fold = fold,
    at_boundary = NA, gap_q1 = NA_real_, gap_q2 = NA_real_,
    stringsAsFactors = FALSE
  )

  fold_results <- silk_parallel_lapply(seq_len(number_folds), function(fold_index) {
    fit_ids <- ids[fold != fold_index]
    hold_ids <- ids[fold == fold_index]
    fit_visits <- train_visits[train_visits$id %in% fit_ids, , drop = FALSE]
    fit_subjects <- train_subjects[train_subjects$id %in% fit_ids, , drop = FALSE]
    hold_visits <- train_visits[train_visits$id %in% hold_ids, , drop = FALSE]
    fit <- fit_registration_multistart(
      fit_visits, fit_subjects, grid,
      seed = if (!is.null(seed)) seed + fold_index else NULL,
      biomarker_kernel = biomarker_kernel
    )
    list(
      prediction = predict_registration_shift(fit$template, hold_visits, grid),
      multistart = fit$multistart,
      fold = fold_index
    )
  })

  fold_multistart <- do.call(rbind, lapply(fold_results, function(res) {
    m <- res$multistart
    if (is.null(m)) return(NULL)
    data.frame(
      fold = res$fold,
      n_starts = m$n_starts,
      n_successful_starts = m$n_successful_starts,
      best_start = m$best_start,
      best_objective = m$best_objective,
      objective_spread = m$objective_spread,
      second_best_gap = m$second_best_gap,
      stringsAsFactors = FALSE
    )
  }))

  for (res in fold_results) {
    prediction <- res$prediction
    row <- match(prediction$id, output$id)
    output$e_hat[row] <- prediction$e_hat
    output$e_clock[row] <- prediction$e_clock
    output$at_boundary[row] <- prediction$at_boundary
    output$gap_q1[row] <- prediction$gap_q1
    output$gap_q2[row] <- prediction$gap_q2
  }
  if (anyNA(output$e_hat)) stop("Cross-fitted registration produced missing shifts.", call. = FALSE)

  subject_order <- match(output$id, train_subjects$id)
  output$S_hat <- train_subjects$A_obs[subject_order] - output$e_hat
  final_fit <- fit_registration_multistart(
    train_visits, train_subjects, grid,
    seed = if (!is.null(seed)) seed + 777L else NULL,
    biomarker_kernel = biomarker_kernel
  )
  structure(
    list(
      train_stage = output,
      final_fit = final_fit,
      final_template = final_fit$template,
      grid = grid,
      biomarker_kernel = final_fit$template$biomarker_builder$kernel,
      characteristic = TRUE,
      implementation = final_fit$implementation,
      # Optimization-stability diagnostics promised by the algorithm section.
      multistart = final_fit$multistart,
      fold_multistart = fold_multistart
    ),
    class = "silk_registration"
  )
}
