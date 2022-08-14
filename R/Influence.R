# Purpose: Check calculation of influence funciton components.
# Updated: 2022-08-13

#' Calculate Mu
#' 
#' Evaluate \eqn{\mu(t; \tau) = \int_{t}^{\tau}{S(u)d(u)/y(u)}du}.
#'
#' @param d Value of d(t) at each time point.
#' @param surv Value of S(t) at each time point.
#' @param unique_times Unique values of time t.
#' @param y Value of y(t) at each time point.
#' @return Numeric vector of \mu(t; tau).
#' @noRd
CalcMu <- function(d, surv, unique_times, y) {
  delta_t <- c(0, diff(unique_times))
  out <- rev(cumsum(rev(delta_t * surv * d / y)))
  return(out)
}


#' Calculate I1
#' 
#' Calculate \eqn{I_{1,i} = \int_{0}^{\tau} -\mu(t; \tau) dM_{i}(t) / y(t)}.
#' 
#' @param dm Matrix of dM_{i}(t).
#' @param mu Vector of \mu(t; tau).
#' @param y Vector of y(t).
#' @return Vector with I1 for each subject.
#' @noRd
CalcI1 <- function(dm, mu, y) {
  n <- nrow(dm)
  out <- lapply(seq_len(n), function(i){
    i1 <- -1 * sum( dm[i, ] * mu / y )
  }) 
  out <- do.call(c, out)
  return(out)
}
