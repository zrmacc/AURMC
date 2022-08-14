# Purpose: Check calculation of influence funciton components.
# Updated: 2022-08-14


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


#' Calculate I2
#' 
#' Calculate \eqn{I_{2,i} = S(t)\{D_{i}(t) - d(t)\}/y(t)}.
#' 
#' @param d Vector of d(t).
#' @param surv Vector of S(t).
#' @param unique_times Vector of unique times t.
#' @param value_mat Matrix of D_{i}(t).
#' @param y Vector of y(t).
#' @return Vector with I2 for each subject.
#' @noRd
CalcI2 <- function(d, surv, unique_times, value_mat, y) {
  n <- nrow(value_mat)
  dt <- c(diff(unique_times), 0)
  out <- lapply(seq_len(n), function(i){
    di <- value_mat[i, ]
    i2 <- sum(surv * (di - d) * dt / y)
  }) 
  out <- do.call(c, out)
  return(out)
}


#' Calculate I3
#' 
#' Calculate \eqn{I_{3,i} = -S(t)d(t)\{Y_{i}(t) - y(t)\} / y^2(t)}.
#'
#' @param d Vector of d(t).
#' @param risk_mat Matrix of Y_{i}(t).
#' @param surv Vector of S(t).
#' @param unique_times Vector of unique times t.
#' @param y Vector of y(t).
#' @return Vector with I3 for each subject.
#' @noRd
CalcI3 <- function(d, risk_mat, surv, unique_times, y) {
  n <- nrow(risk_mat)
  dt <- c(diff(unique_times), 0)
  y2 <- y * y
  out <- lapply(seq_len(n), function(i){
    yi <- risk_mat[i, ]
    i3 <- -1 * sum(surv * d * (yi - y) * dt / y2)
  }) 
  out <- do.call(c, out)
  return(out)
}
