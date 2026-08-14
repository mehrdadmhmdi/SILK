testthat::test_that("the method roster is exact and contains no legacy variants", {
  expected <- c(
    "SILK-Gaussian",
    "SILK-Laplace",
    "SILK-Matern32",
    "Recorded-Cox",
    "Recorded-Beran",
    "MMLM",
    "JM",
    "RSF",
    "DeepSurv",
    "TimeError-Integrated-Landmark",
    "Oracle-Cox",
    "Oracle-Beran"
  )

  defaults <- silk_default_options()
  testthat::expect_identical(defaults$METHODS, expected)
  testthat::expect_identical(defaults$METHOD_ORDER, expected)
  testthat::expect_identical(names(defaults$METHOD_LABELS), expected)
})

testthat::test_that("the four age-only benchmark APIs do not accept visits", {
  testthat::expect_identical(names(formals(fit_recorded_cox)), "train_subjects")
  testthat::expect_identical(names(formals(fit_recorded_beran)), "train_subjects")
  testthat::expect_identical(names(formals(fit_oracle_cox)), "train_subjects")
  testthat::expect_identical(names(formals(fit_oracle_beran)), "train_subjects")

  testthat::expect_false("test_visits" %in% names(formals(predict_recorded_cox)))
  testthat::expect_false("test_visits" %in% names(formals(predict_recorded_beran)))
  testthat::expect_false("test_visits" %in% names(formals(predict_oracle_cox)))
  testthat::expect_false("test_visits" %in% names(formals(predict_oracle_beran)))
})

testthat::test_that("the four benchmarks fit age alone", {
  subjects <- data.frame(
    id = seq_len(8L),
    A_obs = c(42, 45, 49, 52, 56, 59, 63, 67),
    A_star = c(40, 44, 47, 51, 55, 58, 62, 66),
    U = c(2.0, 4.0, 1.5, 3.0, 2.5, 4.5, 3.5, 5.0),
    delta = c(1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L)
  )

  recorded_cox <- fit_recorded_cox(subjects)
  oracle_cox <- fit_oracle_cox(subjects)
  recorded_beran <- fit_recorded_beran(subjects)
  oracle_beran <- fit_oracle_beran(subjects)

  testthat::expect_identical(recorded_cox$fit$type, "cox_age_null")
  testthat::expect_identical(oracle_cox$fit$type, "cox_age_null")
  testthat::expect_length(recorded_cox$fit$keep, 0L)
  testthat::expect_length(oracle_cox$fit$keep, 0L)
  testthat::expect_equal(ncol(recorded_beran$fit$state), 1L)
  testthat::expect_equal(ncol(oracle_beran$fit$state), 1L)
})

testthat::test_that("SILK defaults to the Gaussian method ID", {
  testthat::expect_identical(eval(formals(fit_silk)$method), "SILK-Gaussian")
})
