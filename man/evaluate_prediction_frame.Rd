\name{evaluate_prediction_frame}
\alias{evaluate_prediction_frame}
\title{Evaluate a Complete Prediction Frame}
\description{
  Computes per-horizon metrics, replicate-level metrics, method summaries,
  paired differences, and calibration bins from a prediction frame.
}
\usage{
evaluate_prediction_frame(pred)
}
\arguments{
  \item{pred}{A prediction frame data frame as returned by \code{\link{predict_silk}}.}
}
\value{
  A list with elements: per_horizon, replicate_metrics, summary, paired,
  calibration_bins.
}
\examples{
\dontrun{
dat <- generate_dataset_fixed(200, "mean_moderate", seed = 1)
fit <- fit_silk(dat$subjects, dat$visits, shift_range = c(-12, 12), seed = 1)
test <- generate_dataset_fixed(100, "mean_moderate", seed = 2)
pred <- predict_silk(fit, test$subjects, test$visits)
eval_results <- evaluate_prediction_frame(pred)
}
}
