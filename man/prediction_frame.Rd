\name{prediction_frame}
\alias{prediction_frame}
\title{Build a Prediction Frame}
\usage{prediction_frame(subjects, landmark, horizons, risk_mat, context)}
\arguments{
  \item{subjects}{Data frame with id and D_star columns.}
  \item{landmark}{Numeric vector of landmark ages.}
  \item{horizons}{Numeric vector of prediction horizons.}
  \item{risk_mat}{Matrix of predicted risks (subjects x horizons).}
  \item{context}{Method context list.}
}
\value{Standardized prediction data frame.}
\description{Creates a standardized data frame of risk predictions across subjects and horizons.}
