comparison_env <- new.env(parent = asNamespace("SILK"))
suppressWarnings(
  sys.source(
    testthat::test_path("..", "..", "example", "sim codes", "R", "methods_common.R"),
    envir = comparison_env
  )
)
suppressWarnings(
  sys.source(
    testthat::test_path("..", "..", "example", "sim codes", "R", "methods_mmlm_recorded.R"),
    envir = comparison_env
  )
)
suppressWarnings(
  sys.source(
    testthat::test_path("..", "..", "example", "sim codes", "R", "methods_jm_recorded.R"),
    envir = comparison_env
  )
)

testthat::test_that("misspecified marker differences remove subject-level offsets", {
  visits <- data.frame(
    id = rep(1:3, each = 4),
    visit = rep(1:4, times = 3),
    A_obs_il = rep(24:27, times = 3),
    B1 = c(1, 2, 4, 7, 11, 12, 14, 17, -6, -5, -3, 0)
  )
  shifted <- visits
  shifted$B1 <- shifted$B1 + rep(c(100, -50, 25), each = 4)

  differenced <- comparison_env$marker_values_for_specification(
    visits, "B1", "misspecified"
  )
  differenced_shifted <- comparison_env$marker_values_for_specification(
    shifted, "B1", "misspecified"
  )

  testthat::expect_equal(differenced, differenced_shifted)
  testthat::expect_equal(
    differenced,
    rep(c(NA, 1, 2, 3), 3),
    tolerance = 1e-12
  )
  testthat::expect_identical(
    comparison_env$MISSPECIFIED_ASSOCIATION_SHRINKAGE,
    50
  )
  testthat::expect_identical(
    comparison_env$marker_values_for_specification(visits, "B1", "correct"),
    as.numeric(visits$B1)
  )
})

testthat::test_that("method labels remain unchanged", {
  expected <- c("MMLM-Misspecified", "JM-Misspecified")
  testthat::expect_true(all(expected %in% silk_default_options()$METHODS))
  testthat::expect_identical(
    silk_default_options()$METHODS[match(expected, silk_default_options()$METHODS)],
    expected
  )
})

testthat::test_that("misspecified MMLM uses a finite ridge event model", {
  set.seed(108)
  n <- 120L
  subjects <- data.frame(
    U = stats::rexp(n, rate = 0.2),
    delta = stats::rbinom(n, 1, 0.7)
  )
  x <- cbind(
    X1 = stats::rnorm(n),
    X2 = stats::rbinom(n, 1, 0.45),
    current_biomarker = stats::rnorm(n)
  )
  fit <- comparison_env$fit_mmlm_misspecified_event_model(subjects, x)
  risk <- comparison_env$predict_mmlm_misspecified_event_risk(
    fit, x[1:8, , drop = FALSE], c(1, 2)
  )

  testthat::expect_identical(fit$type, "cox_residual_ridge")
  testthat::expect_identical(fit$theta, 50)
  testthat::expect_equal(dim(risk), c(8L, 2L))
  testthat::expect_true(all(is.finite(risk) & risk > 0 & risk < 1))
})

testthat::test_that("misspecified JM has supported random slopes and prior shape", {
  testthat::skip_if_not_installed("nlme")
  set.seed(109)
  n <- 30L
  visits <- data.frame(
    id = rep(seq_len(n), each = 4L),
    A_obs_il = rep(24:27, times = n),
    X1 = rep(stats::rnorm(n), each = 4L),
    X2 = rep(stats::rbinom(n, 1, 0.45), each = 4L)
  )
  subject_offset <- rep(stats::rnorm(n, sd = 2), each = 4L)
  subject_slope <- rep(stats::rnorm(n, sd = 0.3), each = 4L)
  visits$B1 <- subject_offset + subject_slope * visits$A_obs_il +
    stats::rnorm(nrow(visits), sd = 0.4)

  fit <- comparison_env$fit_jm_longitudinal_model(visits, "misspecified")
  priors <- comparison_env$jm_association_priors("misspecified")

  testthat::expect_s3_class(fit$fit, "lme")
  testthat::expect_match(fit$random_structure, "Gaussian random intercept")
  testthat::expect_s3_class(fit$fit$modelStruct$reStruct[[1]], "pdDiag")
  testthat::expect_identical(fit$marker_transform, "first_difference")
  testthat::expect_true(is.finite(fit$random_time_center))
  testthat::expect_true("A_obs_centered" %in% names(nlme::getData(fit$fit)))
  random_sd <- sqrt(diag(as.matrix(nlme::getVarCov(
    fit$fit, type = "random.effects"
  ))))
  testthat::expect_gte(
    min(random_sd),
    fit$random_effect_sd_floor_ratio * fit$fit$sigma - 1e-10
  )
  testthat::expect_equal(
    mean(nlme::getData(fit$fit)$A_obs_centered), 0,
    tolerance = 1e-12
  )
  testthat::expect_length(priors$mean_alphas, 1L)
  testthat::expect_length(priors$Tau_alphas, 1L)
  testthat::expect_equal(priors$Tau_alphas[[1]], matrix(50, 1L, 1L))
})
