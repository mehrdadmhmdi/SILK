\name{prediction_per_horizon_metrics}
\alias{prediction_per_horizon_metrics}
\alias{prediction_replicate_metrics}
\alias{prediction_method_summary}
\alias{prediction_paired_differences}
\title{Prediction Evaluation Functions}
\usage{
prediction_per_horizon_metrics(pred)
prediction_replicate_metrics(per_horizon)
prediction_method_summary(replicate_metrics)
prediction_paired_differences(replicate_metrics, reference = "SILK")
}
\arguments{
  \item{pred}{A prediction frame data frame.}
  \item{per_horizon}{Output of \code{prediction_per_horizon_metrics}.}
  \item{replicate_metrics}{Output of \code{prediction_replicate_metrics}.}
  \item{reference}{Character. Reference method name.}
}
\value{Data frames of computed metrics at each level of aggregation.}
\description{Compute prediction performance metrics at per-horizon, replicate, and summary levels.}
