# Purpose: Check calculation of influence funciton components.
# Updated: 2022-08-13

#' Calculate Mu
#' 
#' Evaluate $\mu(t; \tau) = \int_{t}^{\tau}{S(u)d(u)/y(u)}du$.
#'
#' @param d Value of d(t) at each time point.
#' @param surv Value of S(t) at each time point.
#' @param time Values of time t.
#' @param y Value of y(t) at each time point.
#' @return Numeric vector of $\mu(t; tau)$.
#' @noRd
CalcMu <- function(
  d,
  surv,
  time,
  y
) {
  delta_t <- c(0, diff(time))
  out <- rev(cumsum(rev(delta_t * surv * d / y)))
  return(out)
}
