\name{clip_probability}
\alias{clip_probability}
\title{Clip Probabilities}
\usage{clip_probability(x, eps = 1e-06)}
\arguments{
  \item{x}{Numeric vector or matrix.}
  \item{eps}{Small positive number for clipping bounds.}
}
\value{Clipped probabilities preserving dimensions.}
\description{Clips probabilities to the interval (eps, 1-eps).}
