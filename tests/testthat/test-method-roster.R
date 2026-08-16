testthat::test_that("package defaults contain only SILK method identifiers", {
  expected <- c(
    "SILK-Gaussian",
    "SILK-Laplace",
    "SILK-Matern32"
  )

  defaults <- silk_default_options()
  testthat::expect_identical(defaults$METHODS, expected)
  testthat::expect_identical(defaults$METHOD_ORDER, expected)
  testthat::expect_identical(names(defaults$METHOD_LABELS), expected)
})

testthat::test_that("survival benchmark wrappers are not part of SILK", {
  testthat::expect_false("fit_recorded_cox" %in% getNamespaceExports("SILK"))
  testthat::expect_false("fit_oracle_cox" %in% getNamespaceExports("SILK"))
})

testthat::test_that("SILK defaults to the Gaussian method ID", {
  testthat::expect_identical(eval(formals(fit_silk)$method), "SILK-Gaussian")
})

testthat::test_that("binary AUC does not overflow for large samples", {
  y <- rep(c(0L, 1L), each = 50000L)
  score <- seq_along(y)
  testthat::expect_equal(binary_auc(y, score), 1)
})

testthat::test_that("the default package feature map has no simulation covariates", {
  subjects <- data.frame(X1 = 1:3, X2 = c(0, 1, 0), X3 = c(-1, 0, 2), X4 = c(2, -1, 0))
  default <- SILK:::base_covariates(subjects, "correct")
  testthat::expect_equal(dim(default), c(nrow(subjects), 0L))
})
