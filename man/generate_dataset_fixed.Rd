\name{generate_dataset_fixed}
\alias{generate_dataset_fixed}
\title{Generate a Simulated Dataset}
\description{
  Generates training or test data under a specified error scenario and visit
  schedule for the SILK simulation framework.
}
\usage{
generate_dataset_fixed(n, scenario_name, schedule_name = NULL, seed = NULL)
}
\arguments{
  \item{n}{Integer. Number of subjects.}
  \item{scenario_name}{Character string. Name of the error scenario.}
  \item{schedule_name}{Character string. Visit schedule name (default: scenario default).}
  \item{seed}{Integer or NULL. Random seed for reproducibility.}
}
\value{
  A list with elements:
  \describe{
    \item{subjects}{Data frame of subject-level information.}
    \item{visits}{Data frame of visit-level biomarker measurements.}
    \item{scenario}{The scenario name used.}
    \item{schedule}{The schedule name used.}
    \item{scenario_row}{The full scenario definition row.}
  }
}
\examples{
dat <- generate_dataset_fixed(100, "mean_moderate", seed = 42)
head(dat$subjects)
head(dat$visits)
}
