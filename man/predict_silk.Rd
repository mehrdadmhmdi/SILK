\name{predict_silk}
\alias{predict_silk}
\title{Predict Risk from a SILK Model}
\description{
  Generates absolute risk predictions for new subjects using a fitted SILK model.
}
\usage{
predict_silk(fit, test_subjects, test_visits, horizons = NULL,
             fold_id = 1L, replicate_id = 1L,
             n_train_setting = NA, time_grid_setting = NA)
}
\arguments{
  \item{fit}{Fitted SILK model from \code{\link{fit_silk}}.}
  \item{test_subjects}{Data frame of test subjects.}
  \item{test_visits}{Data frame of test visits with biomarker columns.}
  \item{horizons}{Numeric vector of prediction horizons. Defaults to \code{silk_opt("PREDICTION_HORIZONS")}.}
  \item{fold_id}{Integer. Fold identifier for bookkeeping.}
  \item{replicate_id}{Integer. Replicate identifier for bookkeeping.}
  \item{n_train_setting}{Numeric. Training set size setting for bookkeeping.}
  \item{time_grid_setting}{Character. Time grid setting for bookkeeping.}
}
\value{
  A prediction frame data frame with columns: subject_id, landmark, horizon,
  method, risk_pred, event_within_horizon, at_risk, fold_id, replicate_id,
  n_train_setting, time_grid_setting.
}
\examples{
\dontrun{
dat <- generate_dataset_fixed(200, "mean_moderate", seed = 1)
fit <- fit_silk(dat$subjects, dat$visits, shift_range = c(-12, 12), seed = 1)
test <- generate_dataset_fixed(100, "mean_moderate", seed = 2)
pred <- predict_silk(fit, test$subjects, test$visits)
}
}
