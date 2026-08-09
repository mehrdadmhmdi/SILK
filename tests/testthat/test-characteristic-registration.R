local_silk_options <- function(...) {
  values <- list(...)
  option_names <- paste0("silk.", names(values))
  old <- options()[option_names]
  missing <- setdiff(option_names, names(old))
  if (length(missing)) old[missing] <- list(NULL)
  withr::defer(options(old), envir = parent.frame())
  do.call(silk_options, values)
}

test_that("Gaussian biomarker kernel is normalized and positive semidefinite", {
  x <- cbind(a = c(-2, -0.5, 1, 3), b = c(1, -1, 0.5, 2))
  kernel <- SILK:::gaussian_rbf_kernel(x, bandwidth = 1.3)
  expect_equal(kernel, t(kernel), tolerance = 1e-12)
  expect_equal(diag(kernel), rep(1, nrow(x)), tolerance = 1e-12)
  expect_gte(min(eigen(kernel, symmetric = TRUE, only.values = TRUE)$values), -1e-10)
})

test_that("Gaussian MMD distinguishes distributions with matching first two moments", {
  p <- matrix(c(-1, -1, 1, 1), ncol = 1)
  q <- matrix(c(-sqrt(2), 0, 0, sqrt(2)), ncol = 1)
  expect_equal(mean(p), mean(q), tolerance = 1e-12)
  expect_equal(stats::var(p), stats::var(q), tolerance = 1e-12)
  mmd_squared <- mean(SILK:::gaussian_rbf_kernel(p, bandwidth = 0.8)) +
    mean(SILK:::gaussian_rbf_kernel(q, bandwidth = 0.8)) -
    2 * mean(SILK:::gaussian_rbf_kernel(p, q, bandwidth = 0.8))
  expect_gt(mmd_squared, 1e-3)
})

test_that("implemented RKHS loss equals its direct kernel expansion", {
  biomarkers <- matrix(c(-1.2, 0.1, 1.4), ncol = 1)
  query <- matrix(0.35, ncol = 1)
  kernel <- SILK:::gaussian_rbf_kernel(biomarkers, bandwidth = 1.1)
  cross_kernel <- SILK:::gaussian_rbf_kernel(biomarkers, query, bandwidth = 1.1)
  A <- matrix(c(0.2, -0.1, 0.3, 0.4, 0.1, -0.2), nrow = 3, ncol = 2)
  psi <- c(0.7, -0.25)
  beta <- as.numeric(A %*% psi)
  direct <- 1 - 2 * sum(beta * cross_kernel) +
    as.numeric(crossprod(beta, kernel %*% beta))
  summary_loss <- SILK:::rkhs_loss_from_summaries(
    matrix(psi, nrow = 1),
    t(crossprod(A, cross_kernel)),
    crossprod(A, kernel %*% A)
  )
  expect_equal(as.numeric(summary_loss), direct, tolerance = 1e-12)
})

test_that("primal template coefficients equal the dual ridge solution", {
  set.seed(4)
  psi_train <- matrix(rnorm(30), nrow = 6, ncol = 5)
  psi_query <- rnorm(5)
  rho <- 0.7
  A <- psi_train %*% solve(crossprod(psi_train) + rho * diag(5))
  primal <- as.numeric(A %*% psi_query)
  input_gram <- tcrossprod(psi_train)
  dual <- as.numeric(solve(input_gram + rho * diag(6), psi_train %*% psi_query))
  expect_equal(primal, dual, tolerance = 1e-10)
})

test_that("missing visits retain nominal schedule positions", {
  subjects <- data.frame(
    id = 1:12,
    A_star = seq(24, 31, length.out = 12),
    eps = 0,
    X1 = 0,
    X2 = 0,
    U_latent = 0
  )
  scenario <- data.frame(
    n_biomarkers = 2L,
    irregular = FALSE,
    missing_rate = 1,
    biomarker_signal = "mean",
    signal_amp = 1,
    sigma_bio = 0.2,
    u_bio_coef = 0
  )
  set.seed(8)
  visits <- SILK:::make_visits(subjects, c(6, 4, 2, 0), scenario)
  by_subject <- split(visits$visit, visits$id)
  expect_true(all(vapply(by_subject, length, integer(1)) == 2L))
  expect_true(all(vapply(by_subject, max, numeric(1)) == 4))
})

test_that("Cox and Beran layers reuse one characteristic registration", {
  local_silk_options(
    N_FOLDS = 2L,
    N_STARTS = 1L,
    N_ALT_ITER_MAX = 1L,
    N_CORES = 1L,
    SHIFT_GRID_STEP = 1,
    TEMPLATE_INPUT_FEATURES = 8L,
    TEMPLATE_QUERY_CHUNK = 100L
  )
  train <- generate_dataset_fixed(30, "no_error", schedule_name = "m2", seed = 31)
  test <- generate_dataset_fixed(12, "no_error", schedule_name = "m2", seed = 32)
  registration <- fit_silk_registration(
    train$subjects, train$visits, shift_grid = c(-1, 0, 1), seed = 33
  )
  expect_s3_class(registration, "silk_registration")
  expect_true(registration$characteristic)
  expect_identical(registration$biomarker_kernel, "gaussian_rbf")

  cox <- fit_silk(train$subjects, train$visits, registration = registration)
  beran <- fit_beran_silk(train$subjects, train$visits, registration = registration)
  expect_identical(cox$registration, registration)
  expect_identical(beran$registration, registration)
  expect_s3_class(predict_silk_registration(registration, test$visits), "data.frame")
  expect_s3_class(predict_silk(cox, test$subjects, test$visits), "data.frame")
  expect_s3_class(predict_beran_silk(beran, test$subjects, test$visits), "data.frame")
})
