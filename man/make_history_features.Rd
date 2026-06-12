\name{make_history_features}
\alias{make_history_features}
\title{Build Landmark History Features}
\usage{make_history_features(subjects, visits, max_biomarkers = 4L)}
\arguments{
  \item{subjects}{Data frame of subjects.}
  \item{visits}{Data frame of visits with biomarker columns.}
  \item{max_biomarkers}{Integer. Maximum number of biomarkers for current values.}
}
\value{Data frame with one row per subject containing history features.}
\description{Creates summary features from biomarker visit histories.}
