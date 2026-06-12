\name{fit_residual_cox}
\alias{fit_residual_cox}
\alias{predict_residual_cox_risk}
\title{Residual Cox Model}
\usage{
fit_residual_cox(subjects, x)
predict_residual_cox_risk(model, x_new, horizons)
}
\arguments{
  \item{subjects}{Data frame with U (follow-up time) and delta (event indicator).}
  \item{x}{Numeric matrix of covariates.}
  \item{model}{Fitted model from \code{fit_residual_cox}.}
  \item{x_new}{Numeric matrix of new covariates.}
  \item{horizons}{Numeric vector of prediction horizons.}
}
\value{
  \code{fit_residual_cox} returns a fitted model object.
  \code{predict_residual_cox_risk} returns a matrix of predicted risks.
}
\description{Fit and predict from a residual-time Cox proportional hazards model.}
