# =============================================================================
# utils-survival.R
# Survival helpers (from surv.R)
# =============================================================================

#' Build a uniform grid from 0 to tau
#' @param tau Numeric. Maximum time.
#' @param n Integer. Number of grid points.
#' @return Numeric vector.
#' @keywords internal
build_u_grid <- function(tau = NULL, n = 80L) {
  if (is.null(tau)) tau <- silk_opt("EVAL_HORIZON")
  seq(0, tau, length.out = n)
}

#' Trapezoidal integration
#' @param x Numeric vector of grid points.
#' @param y Numeric vector of function values.
#' @return Numeric scalar.
#' @keywords internal
trapz <- function(x, y) {
  if (length(x) < 2L) return(0)
  sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(y)])) / 2
}

#' Evaluate a step-function survival curve on a grid
#' @param times Numeric vector of event times.
#' @param surv Numeric vector of survival probabilities.
#' @param grid Numeric vector of evaluation times.
#' @return Numeric vector of survival values on the grid.
#' @keywords internal
step_eval_survival <- function(times, surv, grid) {
  times <- as.numeric(times)
  surv <- as.numeric(surv)
  grid <- as.numeric(grid)
  keep <- is.finite(times) & is.finite(surv)
  if (!any(keep)) return(rep(1, length(grid)))
  z <- data.frame(time = c(0, times[keep]), surv = c(1, surv[keep]))
  z <- z[order(z$time), , drop = FALSE]
  z <- z[!duplicated(z$time, fromLast = TRUE), , drop = FALSE]
  out <- rep(1, length(grid))
  valid_grid <- is.finite(grid)
  if (any(valid_grid)) {
    idx <- findInterval(grid[valid_grid], z$time)
    out[valid_grid] <- z$surv[pmax(idx, 1L)]
  }
  out
}

#' Evaluate a step-function cumulative hazard on a grid
#' @param times Numeric vector of event times.
#' @param cumhaz Numeric vector of cumulative hazard values.
#' @param grid Numeric vector of evaluation times.
#' @return Numeric vector of cumulative hazard values on the grid.
#' @keywords internal
step_eval_cumhaz <- function(times, cumhaz, grid) {
  times <- as.numeric(times)
  cumhaz <- as.numeric(cumhaz)
  grid <- as.numeric(grid)
  keep <- is.finite(times) & is.finite(cumhaz)
  if (!any(keep)) return(rep(0, length(grid)))
  z <- data.frame(time = c(0, times[keep]), cumhaz = c(0, cumhaz[keep]))
  z <- z[order(z$time), , drop = FALSE]
  z <- z[!duplicated(z$time, fromLast = TRUE), , drop = FALSE]
  out <- rep(0, length(grid))
  valid_grid <- is.finite(grid)
  if (any(valid_grid)) {
    idx <- findInterval(grid[valid_grid], z$time)
    out[valid_grid] <- z$cumhaz[pmax(idx, 1L)]
  }
  out
}
