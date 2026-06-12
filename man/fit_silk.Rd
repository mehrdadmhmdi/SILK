\name{fit_silk}
\alias{fit_silk}
\title{Fit a SILK Model}
\description{
  Fits the SILK (Shift-Invariant Learned Kernel) model for absolute risk
  prediction under origin-time measurement error. This includes biomarker
  trajectory registration with cross-fitting and a residual Cox model. For
  real data, the error distribution is not assumed known; the user supplies
  only a candidate shift grid or range over which registration is optimized.
}
\usage{
fit_silk(train_subjects, train_visits, shift_grid = NULL, seed = NULL,
         shift_range = NULL, kernel = NULL, kernel_approx = NULL,
         rff_dim = NULL, rff_seed = seed)
}
\arguments{
  \item{train_subjects}{Data frame of training subjects.}
  \item{train_visits}{Data frame of training visits with biomarker columns.}
  \item{shift_grid}{Numeric vector of candidate origin shifts. If omitted,
  the grid is built from \code{shift_range} and \code{SHIFT_GRID_STEP}.}
  \item{seed}{Integer or NULL. Random seed for reproducibility.}
  \item{shift_range}{Numeric vector of length two giving the lower and upper
  candidate shift values. Defaults to \code{DEFAULT_SHIFT_GRID_MIN} and
  \code{DEFAULT_SHIFT_GRID_MAX}.}
  \item{kernel}{Registration kernel. One of \code{"rbf"}, \code{"matern"},
  \code{"polynomial"}, or \code{"linear"}. Defaults to
  \code{silk_opt("REGISTRATION_KERNEL")}.}
  \item{kernel_approx}{Kernel approximation. Use \code{"exact"} or
  \code{"rff"}. RFF is available for RBF and Matern kernels.}
  \item{rff_dim}{Number of random Fourier features when
  \code{kernel_approx = "rff"}.}
  \item{rff_seed}{Integer or NULL. Random seed for RFF basis generation.}
}
\value{
  A fitted SILK model object (list) for use with \code{\link{predict_silk}}.
}
\examples{
\dontrun{
dat <- generate_dataset_fixed(200, "mean_moderate", seed = 1)
fit <- fit_silk(dat$subjects, dat$visits, shift_range = c(-12, 12), seed = 1)
fit_fast <- fit_silk(dat$subjects, dat$visits, shift_range = c(-12, 12),
                     kernel = "rbf", kernel_approx = "rff",
                     rff_dim = 512, seed = 1)
}
}
