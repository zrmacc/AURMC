# Purpose: Check calculation of influence function components.
# Updated: 2022-08-14


#' Integration Weights
#'
#' @param unique_times Strictly increasing integration grid.
#' @param int_method Integration method.
#' @return Numeric vector of quadrature weights.
#' @noRd
IntegrationWeights <- function(unique_times, int_method = "trapezoid") {
  int_method <- match.arg(int_method, c("left", "right", "trapezoid"))
  delta_t <- diff(unique_times)
  if (length(delta_t) == 0 || any(delta_t <= 0)) {
    stop("At least two strictly increasing time points are required.", call. = FALSE)
  }
  switch(
    int_method,
    left = c(delta_t, 0),
    right = c(0, delta_t),
    trapezoid = {
      weights <- numeric(length(unique_times))
      weights[-length(weights)] <- weights[-length(weights)] + delta_t / 2
      weights[-1] <- weights[-1] + delta_t / 2
      weights
    }
  )
}


#' Calculate Mu
#' 
#' Evaluate the quadrature tail associated with
#' \eqn{\mu(t; \tau) = \int_{t}^{\tau}{S(u-)d(u)/y(u)}du}.
#' The current grid point is excluded because the hazard jump at \eqn{t}
#' affects \eqn{\hat S(u-)} only for \eqn{u > t}.
#'
#' @param d Value of d(t) at each time point.
#' @param surv Value of \eqn{S(t-)} at each time point.
#' @param unique_times Unique values of time t.
#' @param y Value of y(t) at each time point.
#' @return Numeric vector of \mu(t; tau).
#' @noRd
CalcMu <- function(d, surv, unique_times, y, int_method = "trapezoid") {
  weights <- IntegrationWeights(unique_times, int_method)
  inclusive_tail <- rev(cumsum(rev(weights * surv * d / y)))
  out <- c(inclusive_tail[-1], 0)
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
#' Calculate \eqn{I_{2,i} = \int_{0}^{\tau} S(t-)\{D_{i}(t) - d(t)\}/y(t) dt}.
#' 
#' @param d Vector of d(t).
#' @param surv Vector of \eqn{S(t-)}.
#' @param unique_times Vector of unique times t.
#' @param value_mat Matrix of D_{i}(t).
#' @param y Vector of y(t).
#' @return Vector with I2 for each subject.
#' @noRd
CalcI2 <- function(d, surv, unique_times, value_mat, y, int_method = "trapezoid") {
  n <- nrow(value_mat)
  weights <- IntegrationWeights(unique_times, int_method)
  out <- lapply(seq_len(n), function(i){
    di <- value_mat[i, ]
    i2 <- sum(surv * (di - d) * weights / y)
  }) 
  out <- do.call(c, out)
  return(out)
}


#' Calculate I3
#' 
#' Calculate \eqn{I_{3,i} = \int_{0}^{\tau} -S(t-)d(t)\{Y_{i}(t) - y(t)\} /
#' y^2(t) dt}.
#' 
#' @param d Vector of d(t).
#' @param risk_mat Matrix of Y_{i}(t).
#' @param surv Vector of \eqn{S(t-)}.
#' @param unique_times Vector of unique times t.
#' @param y Vector of y(t).
#' @return Vector with I3 for each subject.
#' @noRd
CalcI3 <- function(d, risk_mat, surv, unique_times, y, int_method = "trapezoid") {
  n <- nrow(risk_mat)
  weights <- IntegrationWeights(unique_times, int_method)
  y2 <- y * y
  out <- lapply(seq_len(n), function(i){
    yi <- risk_mat[i, ]
    i3 <- -1 * sum(surv * d * (yi - y) * weights / y2)
  }) 
  out <- do.call(c, out)
  return(out)
}
