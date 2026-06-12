\name{binary_auc}
\alias{binary_auc}
\title{Binary AUC via Mann-Whitney U}
\usage{binary_auc(y, score)}
\arguments{
  \item{y}{Integer vector of binary outcomes (0/1).}
  \item{score}{Numeric vector of predicted scores.}
}
\value{Numeric scalar AUC, or NA if undefined.}
\description{Computes the area under the ROC curve using the Mann-Whitney U statistic.}
