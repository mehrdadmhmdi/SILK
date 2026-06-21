# =============================================================================
# registration.R
# Core registration engine (from reg.R)
# =============================================================================

#' @keywords internal
make_scaler <- function(B) {
  B <- as.matrix(B)
  center <- colMeans(B)
  scale <- apply(B, 2, stats::sd)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  list(center = center, scale = scale)
}

#' @keywords internal
apply_scaler <- function(B, scaler) {
  B <- as.matrix(B)
  B <- sweep(B, 2, scaler$center, "-")
  B <- sweep(B, 2, scaler$scale, "/")
  B
}

#' @keywords internal
make_path_features <- function(B_std, visits_df) {
  B_std <- as.matrix(B_std)
  ids <- visits_df$id
  B_centered <- B_std
  B_delta_landmark <- B_std

  for (id in unique(ids)) {
    idx <- which(ids == id)
    mu_i <- colMeans(B_std[idx, , drop = FALSE])
    B_centered[idx, ] <- sweep(B_std[idx, , drop = FALSE], 2, mu_i, "-")

    lm_idx <- idx[which.min(abs(visits_df$lag[idx]))]
    lm_val <- B_std[lm_idx, ]
    B_delta_landmark[idx, ] <- sweep(B_std[idx, , drop = FALSE], 2, lm_val, "-")
  }

  colnames(B_centered) <- paste0(colnames(B_std), "_subject_centered")
  colnames(B_delta_landmark) <- paste0(colnames(B_std), "_minus_landmark")
  list(centered = B_centered, delta_landmark = B_delta_landmark)
}

#' @keywords internal
make_raw_features <- function(B_std, visits_df, feature_type = c("mean", "silk")) {
  feature_type <- match.arg(feature_type)
  B_std <- as.matrix(B_std)
  path <- make_path_features(B_std, visits_df)

  mean_feats <- cbind(
    B_std,
    0.75 * path$centered,
    0.75 * path$delta_landmark
  )
  if (feature_type == "mean") return(mean_feats)

  p <- ncol(B_std)
  sq <- B_std^2
  colnames(sq) <- paste0(colnames(B_std), "_sq")
  c_sq <- path$centered^2
  colnames(c_sq) <- paste0(colnames(path$centered), "_sq")
  feats <- cbind(mean_feats, 0.45 * sq, 0.35 * c_sq)

  if (p >= 2L) {
    npair <- min(6L, p - 1L)
    inter <- sapply(seq_len(npair), function(k) B_std[, k] * B_std[, k + 1L])
    if (is.null(dim(inter))) inter <- matrix(inter, ncol = 1L)
    colnames(inter) <- paste0("int", seq_len(ncol(inter)))
    feats <- cbind(feats, 0.30 * inter)
  }

  norm2 <- rowMeans(B_std^2)
  centered_norm2 <- rowMeans(path$centered^2)
  feats <- cbind(feats, norm2 = 0.50 * norm2, centered_norm2 = 0.50 * centered_norm2)
  feats
}

#' @keywords internal
make_feature_builder <- function(visits_df, feature_type = c("mean", "silk")) {
  feature_type <- match.arg(feature_type)
  bcols <- bio_columns(visits_df)
  B <- as.matrix(visits_df[, bcols, drop = FALSE])
  biomarker_scaler <- make_scaler(B)
  B_std <- apply_scaler(B, biomarker_scaler)
  F_raw <- make_raw_features(B_std, visits_df, feature_type)
  feature_scaler <- make_scaler(F_raw)
  list(
    feature_type = feature_type,
    biomarker_cols = bcols,
    biomarker_scaler = biomarker_scaler,
    feature_scaler = feature_scaler
  )
}

#' @keywords internal
apply_feature_builder <- function(visits_df, builder) {
  B <- as.matrix(visits_df[, builder$biomarker_cols, drop = FALSE])
  B_std <- apply_scaler(B, builder$biomarker_scaler)
  F_raw <- make_raw_features(B_std, visits_df, builder$feature_type)
  apply_scaler(F_raw, builder$feature_scaler)
}

#' @keywords internal
registration_kernel_config <- function(kernel = NULL, approximation = NULL,
                                       rff_dim = NULL, rff_seed = NULL) {
  if (is.null(kernel)) kernel <- silk_opt("REGISTRATION_KERNEL")
  if (is.null(approximation)) approximation <- silk_opt("REGISTRATION_KERNEL_APPROX")
  if (is.null(rff_dim)) rff_dim <- silk_opt("KERNEL_RFF_DIM")
  if (is.null(rff_seed)) rff_seed <- silk_opt("KERNEL_RFF_SEED")

  kernel <- tolower(as.character(kernel)[1])
  if (kernel %in% c("gaussian", "squared_exponential", "squared-exponential")) {
    kernel <- "rbf"
  }
  if (!kernel %in% c("rbf", "matern", "polynomial", "linear")) {
    stop("Unsupported REGISTRATION_KERNEL: ", kernel, call. = FALSE)
  }

  approximation <- tolower(as.character(approximation)[1])
  if (approximation %in% c("full", "none")) approximation <- "exact"
  if (approximation %in% c("random_fourier", "random-fourier")) approximation <- "rff"
  if (!approximation %in% c("exact", "rff")) {
    stop("Unsupported REGISTRATION_KERNEL_APPROX: ", approximation, call. = FALSE)
  }
  if (identical(approximation, "rff") && !kernel %in% c("rbf", "matern")) {
    stop("RFF approximation is available for rbf and matern kernels only.", call. = FALSE)
  }

  rff_dim <- as.integer(rff_dim)
  if (!is.finite(rff_dim) || rff_dim < 8L) rff_dim <- 8L

  chunk <- as.integer(silk_opt("KERNEL_CHUNK_SIZE"))
  if (!is.finite(chunk) || chunk < 1L) chunk <- 2500L

  list(
    kernel = kernel,
    approximation = approximation,
    rff_dim = rff_dim,
    rff_seed = rff_seed,
    matern_nu = as.numeric(silk_opt("KERNEL_MATERN_NU")),
    polynomial_degree = as.integer(silk_opt("KERNEL_POLYNOMIAL_DEGREE")),
    polynomial_coef0 = as.numeric(silk_opt("KERNEL_POLYNOMIAL_COEF0")),
    polynomial_scale = as.numeric(silk_opt("KERNEL_POLYNOMIAL_SCALE")),
    polynomial_nonnegative = isTRUE(silk_opt("KERNEL_POLYNOMIAL_NONNEGATIVE")),
    chunk_size = chunk
  )
}

#' @keywords internal
registration_scaled_inputs <- function(age, x1) {
  h_age <- silk_opt("H_A_BANDWIDTH")
  h_x <- silk_opt("H_X_BANDWIDTH")
  if (!is.finite(h_age) || h_age <= 0) h_age <- 1
  if (!is.finite(h_x) || h_x <= 0) h_x <- 1
  cbind(age = as.numeric(age) / h_age, x1 = as.numeric(x1) / h_x)
}

#' @keywords internal
pairwise_sqdist <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  out <- outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y)
  pmax(out, 0)
}

#' @keywords internal
matern_from_distance <- function(r, nu) {
  if (!is.finite(nu) || nu <= 0) nu <- 1.5
  if (abs(nu - 0.5) < 1e-8) return(exp(-r))
  if (abs(nu - 1.5) < 1e-8) {
    z <- sqrt(3) * r
    return((1 + z) * exp(-z))
  }
  if (abs(nu - 2.5) < 1e-8) {
    z <- sqrt(5) * r
    return((1 + z + z^2 / 3) * exp(-z))
  }
  z <- sqrt(2 * nu) * r
  out <- (2^(1 - nu) / gamma(nu)) * z^nu * besselK(z, nu)
  out[r <= 1e-12] <- 1
  out[!is.finite(out)] <- 0
  out
}

#' @keywords internal
registration_kernel_matrix <- function(train_x, query_x, config) {
  train_x <- as.matrix(train_x)
  query_x <- as.matrix(query_x)

  if (identical(config$kernel, "rbf")) {
    K <- exp(-0.5 * pairwise_sqdist(train_x, query_x))
  } else if (identical(config$kernel, "matern")) {
    K <- matern_from_distance(sqrt(pairwise_sqdist(train_x, query_x)), config$matern_nu)
  } else if (identical(config$kernel, "polynomial")) {
    deg <- config$polynomial_degree
    if (!is.finite(deg) || deg < 1L) deg <- 2L
    base <- config$polynomial_coef0 + config$polynomial_scale * tcrossprod(train_x, query_x)
    if (isTRUE(config$polynomial_nonnegative)) base <- pmax(base, 0)
    K <- base^deg
  } else {
    K <- config$polynomial_coef0 + config$polynomial_scale * tcrossprod(train_x, query_x)
    K <- pmax(K, 0)
  }

  K[!is.finite(K)] <- 0
  K
}

#' @keywords internal
with_preserved_seed <- function(seed, expr) {
  if (is.null(seed) || length(seed) == 0L || !is.finite(seed)) return(force(expr))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(expr)
}

#' @keywords internal
make_rff_model <- function(input_dim, config, stream = 0L) {
  seed <- suppressWarnings(as.numeric(config$rff_seed))
  if (length(seed) == 0L || !is.finite(seed)) seed <- NULL
  if (!is.null(seed)) seed <- seed + as.integer(stream)

  with_preserved_seed(seed, {
    D <- as.integer(config$rff_dim)
    if (identical(config$kernel, "matern")) {
      nu <- config$matern_nu
      if (!is.finite(nu) || nu <= 0) nu <- 1.5
      df <- 2 * nu
      z <- matrix(stats::rnorm(input_dim * D), nrow = input_dim, ncol = D)
      chi <- pmax(stats::rchisq(D, df = df), 1e-12)
      omega <- sweep(z, 2, sqrt(df / chi), "*")
    } else {
      omega <- matrix(stats::rnorm(input_dim * D), nrow = input_dim, ncol = D)
    }
    list(
      omega = omega,
      phase = stats::runif(D, min = 0, max = 2 * pi),
      n_features = D
    )
  })
}

#' @keywords internal
rff_features <- function(x, model) {
  x <- as.matrix(x)
  z <- x %*% model$omega
  z <- sweep(z, 2, model$phase, "+")
  sqrt(2 / model$n_features) * cos(z)
}

#' @keywords internal
make_template_entry <- function(age, x1, Phi, kernel_config, stream = 0L) {
  input <- registration_scaled_inputs(age, x1)
  entry <- list(
    age = age,
    X1 = x1,
    Phi = Phi,
    input = input,
    kernel = kernel_config
  )
  if (identical(kernel_config$approximation, "rff")) {
    model <- make_rff_model(ncol(input), kernel_config, stream)
    Z <- rff_features(input, model)
    entry$rff_model <- model
    entry$rff_phi <- crossprod(Z, Phi)
    entry$rff_weight <- colSums(Z)
  }
  entry
}

#' @keywords internal
anchor_basis <- function(subjects, ids, mode = NULL) {
  if (is.null(mode)) mode <- silk_opt("ANCHOR_MODE")
  ss <- subjects[match(ids, subjects$id), , drop = FALSE]
  if (mode == "intercept") {
    B <- matrix(1, nrow = length(ids), ncol = 1L)
    colnames(B) <- "intercept"
    return(B)
  }
  x <- ss$X1
  x <- as.numeric(scale(x))
  x[!is.finite(x)] <- 0
  B <- cbind(intercept = 1, X1 = x)
  B
}

#' @keywords internal
project_to_anchor <- function(e, A, lower, upper) {
  A <- as.matrix(A)
  if (ncol(A) == 0L) return(pmin(pmax(e, lower), upper))
  G <- crossprod(A)
  rhs <- colSums(A * e)
  adj <- tryCatch(A %*% solve(G + 1e-8 * diag(ncol(A)), rhs), error = function(err) 0)
  e2 <- as.numeric(e - adj)
  pmin(pmax(e2, lower), upper)
}

#' @keywords internal
constrained_grid_update <- function(L_mat, grid, anchor_mat) {
  n <- nrow(L_mat)
  A <- as.matrix(anchor_mat)
  if (nrow(A) != n) stop("anchor_mat has wrong number of rows", call. = FALSE)
  lower <- min(grid); upper <- max(grid)
  SHIFT_GRID_STEP <- silk_opt("SHIFT_GRID_STEP")

  choose_for_lambda <- function(lambda) {
    score <- as.vector(A %*% lambda)
    adj <- L_mat + outer(score, grid, "*")
    grid[max.col(-adj, ties.method = "first")]
  }

  if (ncol(A) == 1L) {
    fsum <- function(lambda) sum(A[, 1] * choose_for_lambda(lambda))
    lo <- -1; hi <- 1
    flo <- fsum(lo); fhi <- fsum(hi)
    expand <- 0L
    while (is.finite(flo) && is.finite(fhi) && flo * fhi > 0 && expand < 40L) {
      lo <- lo * 2; hi <- hi * 2
      flo <- fsum(lo); fhi <- fsum(hi)
      expand <- expand + 1L
    }
    if (is.finite(flo) && is.finite(fhi) && flo * fhi <= 0) {
      for (ii in seq_len(60L)) {
        mid <- (lo + hi) / 2
        fm <- fsum(mid)
        if (abs(fm) <= max(1e-8, 0.01 * n * SHIFT_GRID_STEP)) break
        if (flo * fm <= 0) {
          hi <- mid; fhi <- fm
        } else {
          lo <- mid; flo <- fm
        }
      }
      e <- choose_for_lambda((lo + hi) / 2)
    } else {
      e <- grid[apply(L_mat, 1L, which.min)]
    }
    return(project_to_anchor(e, A, lower, upper))
  }

  obj_lambda <- function(lambda) {
    e <- choose_for_lambda(lambda)
    m <- colSums(A * e)
    sum(m^2) + 1e-8 * sum(lambda^2)
  }

  starts <- list(rep(0, ncol(A)))
  for (s in seq_len(4L)) starts[[length(starts) + 1L]] <- stats::rnorm(ncol(A), sd = 0.5 * s)

  best <- NULL; bestval <- Inf
  for (st in starts) {
    fit <- tryCatch(stats::optim(st, obj_lambda, method = "Nelder-Mead", control = list(maxit = 200)),
                    error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value) && fit$value < bestval) {
      best <- fit$par; bestval <- fit$value
    }
  }
  if (is.null(best)) best <- rep(0, ncol(A))
  e <- choose_for_lambda(best)
  project_to_anchor(e, A, lower, upper)
}

#' @keywords internal
build_position_template <- function(visits_df, subjects, e_vec, ids, feature_type,
                                    kernel_config = NULL) {
  if (is.null(kernel_config)) kernel_config <- registration_kernel_config()
  builder <- make_feature_builder(visits_df, feature_type)
  Phi <- apply_feature_builder(visits_df, builder)
  id_index <- match(visits_df$id, ids)
  cal_age <- visits_df$A_obs_il - e_vec[id_index]
  visits_df$.cal_age <- cal_age
  visits_df$.row <- seq_len(nrow(visits_df))

  templates <- list()
  visit_stream <- 0L
  for (vv in sort(unique(visits_df$visit))) {
    visit_stream <- visit_stream + 1L
    idx <- which(visits_df$visit == vv)
    templates[[as.character(vv)]] <- c(
      list(visit = vv),
      make_template_entry(
        visits_df$.cal_age[idx],
        visits_df$X1[idx],
        Phi[idx, , drop = FALSE],
        kernel_config,
        stream = visit_stream
      )
    )
  }

  list(
    builder = builder,
    templates = templates,
    ids = ids,
    feature_type = feature_type,
    kernel = kernel_config
  )
}

#' @keywords internal
predict_template_mu <- function(tmp_v, q_age, q_x1) {
  kernel_config <- tmp_v$kernel
  if (is.null(kernel_config)) kernel_config <- registration_kernel_config()
  M <- length(q_age)
  d <- ncol(tmp_v$Phi)
  out <- matrix(NA_real_, nrow = M, ncol = d)
  chunk <- kernel_config$chunk_size
  for (st in seq(1L, M, by = chunk)) {
    en <- min(M, st + chunk - 1L)
    qa <- q_age[st:en]
    qx <- q_x1[st:en]
    query_x <- registration_scaled_inputs(qa, qx)

    if (identical(kernel_config$approximation, "rff")) {
      Zq <- rff_features(query_x, tmp_v$rff_model)
      den <- as.numeric(Zq %*% tmp_v$rff_weight)
      num <- Zq %*% tmp_v$rff_phi
      ok <- is.finite(den) & den > 1e-10
      block <- matrix(NA_real_, nrow = length(qa), ncol = d)
      if (any(ok)) block[ok, ] <- sweep(num[ok, , drop = FALSE], 1, den[ok], "/")
      if (any(!ok)) {
        W <- registration_kernel_matrix(tmp_v$input, query_x[!ok, , drop = FALSE], kernel_config)
        rs <- colSums(W)
        rs[!is.finite(rs) | rs <= 1e-12] <- 1e-12
        block[!ok, ] <- sweep(t(W) %*% tmp_v$Phi, 1, rs, "/")
      }
      out[st:en, ] <- block
    } else {
      W <- registration_kernel_matrix(tmp_v$input, query_x, kernel_config)
      rs <- colSums(W)
      rs[!is.finite(rs) | rs <= 1e-12] <- 1e-12
      out[st:en, ] <- sweep(t(W) %*% tmp_v$Phi, 1, rs, "/")
    }
  }
  out
}

#' @keywords internal
loss_grid_from_template <- function(template, visits_df, grid) {
  ids <- sort(unique(visits_df$id))
  n <- length(ids)
  G <- length(grid)
  Phi <- apply_feature_builder(visits_df, template$builder)
  visits_df$.phi_row <- seq_len(nrow(visits_df))

  L_total <- matrix(0, nrow = n, ncol = G)
  N_count <- matrix(0, nrow = n, ncol = G)

  for (vv in sort(unique(visits_df$visit))) {
    idx <- which(visits_df$visit == vv)
    if (!length(idx)) next
    tmp_v <- template$templates[[as.character(vv)]]
    if (is.null(tmp_v)) next

    sv <- visits_df[idx, , drop = FALSE]
    nr <- nrow(sv)
    q_age <- as.vector(outer(sv$A_obs_il, grid, "-"))
    q_x1 <- rep(sv$X1, times = G)
    mu <- predict_template_mu(tmp_v, q_age, q_x1)
    phi_rep <- Phi[idx, , drop = FALSE][rep(seq_len(nr), times = G), , drop = FALSE]
    sq <- rowMeans((phi_rep - mu)^2)
    sq_mat <- matrix(sq, nrow = nr, ncol = G)
    row_id <- match(sv$id, ids)
    for (gg in seq_len(G)) {
      rs_loss <- rowsum(sq_mat[, gg], row_id, reorder = FALSE)
      rs_n <- rowsum(rep(1, nr), row_id, reorder = FALSE)
      acc_loss <- numeric(n)
      acc_n <- numeric(n)
      acc_loss[as.integer(rownames(rs_loss))] <- as.vector(rs_loss)
      acc_n[as.integer(rownames(rs_n))] <- as.vector(rs_n)
      L_total[, gg] <- L_total[, gg] + acc_loss
      N_count[, gg] <- N_count[, gg] + acc_n
    }
  }

  N_count[N_count <= 0] <- 1
  L_total / N_count
}

#' @keywords internal
registration_total_objective <- function(template, visits_df, e_vec, grid_dummy = NULL) {
  ids <- sort(unique(visits_df$id))
  Phi <- apply_feature_builder(visits_df, template$builder)
  loss <- numeric(length(ids)); count <- numeric(length(ids))
  for (vv in sort(unique(visits_df$visit))) {
    idx <- which(visits_df$visit == vv)
    tmp_v <- template$templates[[as.character(vv)]]
    if (is.null(tmp_v)) next
    sv <- visits_df[idx, , drop = FALSE]
    row_id <- match(sv$id, ids)
    e <- e_vec[row_id]
    mu <- predict_template_mu(tmp_v, sv$A_obs_il - e, sv$X1)
    sq <- rowMeans((Phi[idx, , drop = FALSE] - mu)^2)
    rs_loss <- rowsum(sq, row_id, reorder = FALSE)
    rs_n <- rowsum(rep(1, length(sq)), row_id, reorder = FALSE)
    acc_loss <- numeric(length(ids))
    acc_n <- numeric(length(ids))
    acc_loss[as.integer(rownames(rs_loss))] <- as.vector(rs_loss)
    acc_n[as.integer(rownames(rs_n))] <- as.vector(rs_n)
    loss <- loss + acc_loss
    count <- count + acc_n
  }
  mean(loss / pmax(count, 1))
}

#' @keywords internal
fit_registration_train <- function(visits_df, subjects, feature_type = c("mean", "silk"),
                                   init_e = NULL, grid = NULL,
                                   max_iter = NULL, tol = NULL,
                                   kernel_config = NULL) {
  feature_type <- match.arg(feature_type)
  if (is.null(kernel_config)) kernel_config <- registration_kernel_config()
  if (is.null(grid)) grid <- seq(silk_opt("DEFAULT_SHIFT_GRID_MIN"), silk_opt("DEFAULT_SHIFT_GRID_MAX"), by = silk_opt("SHIFT_GRID_STEP"))
  if (is.null(max_iter)) max_iter <- silk_opt("N_ALT_ITER_MAX")
  if (is.null(tol)) tol <- silk_opt("ALT_TOL")
  SHIFT_GRID_STEP <- silk_opt("SHIFT_GRID_STEP")
  ANCHOR_MODE <- silk_opt("ANCHOR_MODE")

  ids <- sort(unique(visits_df$id))
  n <- length(ids)
  if (is.null(init_e)) init_e <- rep(0, n)
  stopifnot(length(init_e) == n)
  A <- anchor_basis(subjects, ids, ANCHOR_MODE)
  e_vec <- project_to_anchor(as.numeric(init_e), A, min(grid), max(grid))
  obj_trace <- numeric(0)

  for (iter in seq_len(max_iter)) {
    template <- build_position_template(visits_df, subjects, e_vec, ids, feature_type, kernel_config)
    L <- loss_grid_from_template(template, visits_df, grid)
    new_e <- constrained_grid_update(L, grid, A)
    obj <- registration_total_objective(template, visits_df, new_e)
    obj_trace <- c(obj_trace, obj)
    max_change <- max(abs(new_e - e_vec))
    e_vec <- new_e
    if (length(obj_trace) >= 2L && abs(diff(utils::tail(obj_trace, 2))) < tol && max_change < SHIFT_GRID_STEP) break
  }

  template <- build_position_template(visits_df, subjects, e_vec, ids, feature_type, kernel_config)
  list(e_train = e_vec, template = template, ids = ids, obj_trace = obj_trace, grid = grid)
}

#' @keywords internal
fit_registration_multistart <- function(visits_df, subjects, feature_type = c("mean", "silk"),
                                        mean_reg_init = NULL, grid = NULL, seed = NULL,
                                        kernel_config = NULL) {
  feature_type <- match.arg(feature_type)
  if (is.null(kernel_config)) kernel_config <- registration_kernel_config()
  if (is.null(grid)) grid <- seq(silk_opt("DEFAULT_SHIFT_GRID_MIN"), silk_opt("DEFAULT_SHIFT_GRID_MAX"), by = silk_opt("SHIFT_GRID_STEP"))
  N_STARTS <- silk_opt("N_STARTS")
  RANDOM_START_SD <- silk_opt("RANDOM_START_SD")
  ANCHOR_MODE <- silk_opt("ANCHOR_MODE")

  if (!is.null(seed)) set.seed(seed)
  ids <- sort(unique(visits_df$id))
  n <- length(ids)
  A <- anchor_basis(subjects, ids, ANCHOR_MODE)

  starts <- list(zero = rep(0, n))
  if (N_STARTS >= 2L) starts$random <- project_to_anchor(stats::rnorm(n, sd = RANDOM_START_SD), A, min(grid), max(grid))
  if (N_STARTS >= 3L && !is.null(mean_reg_init)) {
    starts$meanreg <- project_to_anchor(mean_reg_init, A, min(grid), max(grid))
  } else if (N_STARTS >= 3L) {
    starts$random2 <- project_to_anchor(stats::rnorm(n, sd = 2 * RANDOM_START_SD), A, min(grid), max(grid))
  }

  best <- NULL; best_obj <- Inf; best_start <- NA_character_
  for (nm in names(starts)) {
    fit <- fit_registration_train(
      visits_df, subjects, feature_type, starts[[nm]], grid,
      kernel_config = kernel_config
    )
    obj <- registration_total_objective(fit$template, visits_df, fit$e_train)
    if (is.finite(obj) && obj < best_obj) {
      best <- fit; best_obj <- obj; best_start <- nm
    }
  }
  best$best_obj <- best_obj
  best$best_start <- best_start
  best
}

#' @keywords internal
refine_profile_min <- function(loss_row, grid) {
  SHIFT_GRID_STEP <- silk_opt("SHIFT_GRID_STEP")
  PROFILE_TEMPERATURE <- silk_opt("PROFILE_TEMPERATURE")
  PROFILE_LOCAL_RADIUS <- silk_opt("PROFILE_LOCAL_RADIUS")

  j <- which.min(loss_row)
  ghat <- grid[j]

  win <- which(abs(grid - ghat) <= max(2 * SHIFT_GRID_STEP, PROFILE_LOCAL_RADIUS))
  if (length(win) >= 3L && length(unique(grid[win])) >= 3L) {
    fit <- tryCatch(stats::lm(loss_row[win] ~ grid[win] + I(grid[win]^2)), error = function(e) NULL)
    if (!is.null(fit)) {
      cf <- stats::coef(fit)
      if (length(cf) == 3L && is.finite(cf[2]) && is.finite(cf[3]) && cf[3] > 0) {
        v <- -cf[2] / (2 * cf[3])
        if (is.finite(v) && v >= min(grid[win]) && v <= max(grid[win])) ghat <- v
      }
    }
  }

  if (is.finite(PROFILE_TEMPERATURE) && PROFILE_TEMPERATURE > 0) {
    loc <- which(abs(grid - grid[j]) <= PROFILE_LOCAL_RADIUS)
    rel <- loss_row[loc] - min(loss_row[loc], na.rm = TRUE)
    if (is.finite(max(rel, na.rm = TRUE)) && max(rel, na.rm = TRUE) > 2 * PROFILE_TEMPERATURE) {
      w <- exp(-rel / PROFILE_TEMPERATURE)
      if (sum(w) > 0 && all(is.finite(w))) ghat <- sum(grid[loc] * w) / sum(w)
    }
  }

  pmin(pmax(ghat, min(grid)), max(grid))
}

#' @keywords internal
predict_registration_shift <- function(template, visits_df, grid) {
  ids <- sort(unique(visits_df$id))
  L <- loss_grid_from_template(template, visits_df, grid)
  best <- max.col(-L, ties.method = "first")
  e_refined <- vapply(seq_along(ids), function(i) refine_profile_min(L[i, ], grid), numeric(1))

  data.frame(
    id = ids,
    e_hat = e_refined,
    e_hat_grid = grid[best],
    min_loss = L[cbind(seq_along(ids), best)],
    at_boundary = best %in% c(1L, length(grid)),
    gap_q1 = vapply(seq_along(ids), function(i) {
      far <- abs(grid - grid[best[i]]) >= 1
      if (!any(far)) return(NA_real_)
      min(L[i, far]) - L[i, best[i]]
    }, numeric(1)),
    gap_q2 = vapply(seq_along(ids), function(i) {
      far <- abs(grid - grid[best[i]]) >= 2
      if (!any(far)) return(NA_real_)
      min(L[i, far]) - L[i, best[i]]
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
crossfit_registration <- function(train_visits, train_subjects, feature_type = c("mean", "silk"),
                                  grid, mean_reg_init_full = NULL, seed = NULL,
                                  kernel_config = NULL) {
  feature_type <- match.arg(feature_type)
  if (is.null(kernel_config)) kernel_config <- registration_kernel_config()
  N_FOLDS <- silk_opt("N_FOLDS")

  if (!is.null(seed)) set.seed(seed)
  ids <- sort(unique(train_visits$id))
  n <- length(ids)
  nfold <- min(N_FOLDS, n)
  fold <- sample(rep(seq_len(nfold), length.out = n))
  out <- data.frame(id = ids, e_hat = NA_real_, S_hat = NA_real_, fold = fold,
                    at_boundary = NA, gap_q1 = NA_real_, gap_q2 = NA_real_,
                    stringsAsFactors = FALSE)

  for (ff in seq_len(nfold)) {
    fit_ids <- ids[fold != ff]
    hold_ids <- ids[fold == ff]
    v_fit <- train_visits[train_visits$id %in% fit_ids, , drop = FALSE]
    s_fit <- train_subjects[train_subjects$id %in% fit_ids, , drop = FALSE]
    v_hold <- train_visits[train_visits$id %in% hold_ids, , drop = FALSE]

    init <- NULL
    if (!is.null(mean_reg_init_full)) init <- mean_reg_init_full[match(sort(unique(v_fit$id)), ids)]

    fit <- fit_registration_multistart(
      v_fit, s_fit, feature_type, init, grid,
      seed = if (!is.null(seed)) seed + ff else NULL,
      kernel_config = kernel_config
    )
    pred <- predict_registration_shift(fit$template, v_hold, grid)
    ii <- match(pred$id, out$id)
    out$e_hat[ii] <- pred$e_hat
    out$at_boundary[ii] <- pred$at_boundary
    out$gap_q1[ii] <- pred$gap_q1
    out$gap_q2[ii] <- pred$gap_q2
  }

  s_order <- match(out$id, train_subjects$id)
  out$S_hat <- train_subjects$A_obs[s_order] - out$e_hat

  final_init <- NULL
  if (!is.null(mean_reg_init_full)) final_init <- mean_reg_init_full[match(ids, ids)]
  final_fit <- fit_registration_multistart(
    train_visits, train_subjects, feature_type, final_init, grid,
    seed = if (!is.null(seed)) seed + 777L else NULL,
    kernel_config = kernel_config
  )

  list(train_stage = out, final_fit = final_fit, final_template = final_fit$template)
}
