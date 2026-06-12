\name{calibration_stats}
\alias{calibration_stats}
\title{Calibration Statistics}
\usage{calibration_stats(y, p)}
\arguments{
  \item{y}{Integer vector of binary outcomes.}
  \item{p}{Numeric vector of predicted probabilities.}
}
\value{Named numeric vector with calibration_intercept, calibration_slope, and calibration_in_large.}
\description{Computes calibration intercept, slope, and calibration-in-the-large.}
