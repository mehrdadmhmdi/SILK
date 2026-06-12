\name{silk_options}
\alias{silk_options}
\alias{silk_opt}
\alias{silk_default_options}
\title{SILK Package Options}
\description{
  Get and set SILK package configuration options.
  \code{silk_opt} retrieves a single option.
  \code{silk_options} sets one or more options.
  \code{silk_default_options} returns the full list of defaults.
}
\details{
  Registration kernel options include \code{REGISTRATION_KERNEL}
  (\code{"rbf"}, \code{"matern"}, \code{"polynomial"}, or \code{"linear"}),
  \code{REGISTRATION_KERNEL_APPROX} (\code{"exact"}, \code{"none"}, or
  \code{"rff"}), \code{KERNEL_RFF_DIM}, \code{KERNEL_RFF_SEED},
  \code{KERNEL_MATERN_NU}, and polynomial-kernel controls. \code{"none"} is
  accepted as an alias for exact calculation. RFF acceleration is available for
  RBF and Matern registration kernels.
}
\usage{
silk_opt(name)
silk_options(...)
silk_default_options()
}
\arguments{
  \item{name}{Character string. The option name.}
  \item{...}{Named arguments of option name-value pairs.}
}
\value{
  \code{silk_opt} returns the current value of the named option.
  \code{silk_options} returns invisible NULL.
  \code{silk_default_options} returns a named list of all default values.
}
\examples{
silk_opt("SHIFT_GRID_STEP")
silk_options(SHIFT_GRID_STEP = 0.1)
silk_opt("SHIFT_GRID_STEP")
# Reset
silk_options(SHIFT_GRID_STEP = 0.2)

# Use random Fourier features for a large real-data registration run
silk_options(REGISTRATION_KERNEL = "rbf",
             REGISTRATION_KERNEL_APPROX = "rff",
             KERNEL_RFF_DIM = 512)

# Use the full exact kernel calculation
silk_options(REGISTRATION_KERNEL = "rbf",
             REGISTRATION_KERNEL_APPROX = "none")
}
