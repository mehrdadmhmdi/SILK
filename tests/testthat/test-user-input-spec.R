testthat::test_that("user data specifications support arbitrary columns and engineered covariates", {
  subjects <- data.frame(
    subject = 1:3,
    age_recorded = c(50, 51, 52),
    follow_up = c(2, 3, 4),
    status = c(1, 0, 1),
    X3 = c(-1, 0, 2),
    stringsAsFactors = FALSE
  )
  visits <- data.frame(
    subject = rep(1:3, each = 2),
    visit_no = rep(1:2, 3),
    age_visit = rep(c(48, 50), 3),
    cd4 = seq_len(6),
    cd8 = rev(seq_len(6)),
    stringsAsFactors = FALSE
  )
  spec <- silk_data_spec(
    subject_id = "subject", recorded_age = "age_recorded",
    event_time = "follow_up", event = "status",
    visit_number = "visit_no", visit_age = "age_visit", visit_lag = NULL,
    biomarker_cols = c("cd4", "cd8"),
    covariate_cols = "X3",
    engineered_covariates = list(X3_sq = function(data) data[["X3"]]^2 - 1),
    template_input_covariates = c("X3_sq", "lag"),
    template_input_bandwidths = c(X3_sq = 1, lag = 1)
  )
  prepared <- prepare_silk_inputs(subjects, visits, spec)
  normalized_subjects <- prepared[["subjects"]]
  normalized_visits <- prepared[["visits"]]
  testthat::expect_equal(normalized_subjects[["A_obs"]], subjects[["age_recorded"]])
  testthat::expect_equal(normalized_subjects[["X3_sq"]], subjects[["X3"]]^2 - 1)
  testthat::expect_equal(normalized_visits[["lag"]], c(2, 0, 3, 1, 4, 2))
  testthat::expect_equal(normalized_visits[["X3_sq"]], rep(subjects[["X3"]]^2 - 1, each = 2))
  testthat::expect_equal(prepared[["spec"]][["anchor_mode"]], "intercept")
})

testthat::test_that("omitted optional inputs stay self-contained", {
  spec <- silk_data_spec(
    subject_id = "subject_id", recorded_age = "age_recorded",
    event_time = "follow_up", event = "status",
    visit_number = "visit_no", visit_age = "age_visit"
  )
  testthat::expect_identical(spec$covariate_cols, character(0))
  testthat::expect_identical(spec$clock_covariates, character(0))
  testthat::expect_identical(spec$engineered_covariates, list())
  testthat::expect_identical(spec$anchor_mode, "intercept")
})
