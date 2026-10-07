# Purpose: Simulate data.
# Updated: 2024-11-09


#' Simulate Subject
#' 
#' @param censoring_rate Rate for the time to censoring.
#' @param death_rate Rate for the time to death.
#' @param idx Subject index.
#' @param tau Truncation time.
#' @param last_missing Should the last value be missing? Default: FALSE.
#' @param value_mean Mean of measurement.
#' @param value_sd Standard deviation of measurement.
#' @return Data.frame.
#' @noRd
SimSubj <- function(
    censoring_rate,
    death_rate,
    idx,
    tau,
    last_missing = FALSE,
    value_mean = 0,
    value_sd = 1
) {
  
  # Censoring time.
  cens <- stats::rexp(n = 1, rate = censoring_rate)
  cens <- min(cens, tau)
  
  # Death time.
  death <- stats::rexp(n = 1, rate = death_rate)
  
  # Status.
  if (death <= cens) {
    status <- 2
  } else {
    status <- 0
  }
  time <- min(death, cens)
  
  # Observations.
  obs_times <- seq(from = 0, to = floor(time))
  obs_values <- stats::rnorm(
    n = length(obs_times) + 1,
    mean = value_mean,
    sd = value_sd
  )
  
  # Output.
  out <- data.frame(
    idx = idx,
    time = c(obs_times, time),
    status = c(rep(1, length(obs_times)), status),
    value = obs_values
  )
  if (last_missing) {
    out$value[nrow(out)] <- NA
  }
  return(out)
}


#' Simulate Repeated Measures Data
#' 
#' Simulates repeated measures data at regular time points subject to
#' exponential censoring and death.
#' 
#' @param censoring_rate Rate for the time to censoring.
#' @param death_rate Rate for the time to death.
#' @param n Number of subjects.
#' @param tau Truncation time.
#' @param last_missing Should the last value be missing? Default: FALSE.
#' @param value_mean Mean of measurement.
#' @param value_sd Standard deviation of measurement.
#' @return A data.frame with subject index, observation time, status, and
#'   repeated measurement value.
#' @examples
#' set.seed(1)
#' GenData(censoring_rate = 0.5, death_rate = 0.25, n = 5, tau = 2)
#' @export
GenData <- function(
    censoring_rate,
    death_rate,
    n,
    tau,
    last_missing = FALSE,
    value_mean = 0,
    value_sd = 1
) {
  rates <- c(censoring_rate, death_rate)
  if (length(rates) != 2 || anyNA(rates) || any(!is.finite(rates)) || any(rates < 0)) {
    stop("Event rates must be finite, non-negative numbers.", call. = FALSE)
  }
  if (length(n) != 1 || !is.finite(n) || n < 1 || n != as.integer(n)) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  if (length(tau) != 1 || !is.finite(tau) || tau <= 0) {
    stop("`tau` must be a single positive number.", call. = FALSE)
  }
  if (length(value_sd) != 1 || !is.finite(value_sd) || value_sd < 0) {
    stop("`value_sd` must be a single non-negative number.", call. = FALSE)
  }
  if (length(value_mean) != 1 || !is.finite(value_mean)) {
    stop("`value_mean` must be a single finite number.", call. = FALSE)
  }
  if (length(last_missing) != 1 || !is.logical(last_missing) || is.na(last_missing)) {
    stop("`last_missing` must be TRUE or FALSE.", call. = FALSE)
  }
  
  data <- lapply(seq_len(n), function(i) {
    subj <- SimSubj(
      censoring_rate = censoring_rate,
      death_rate = death_rate,
      idx = i, 
      tau = tau, 
      last_missing = last_missing, 
      value_mean = value_mean, 
      value_sd = value_sd
    )
    return(subj)
    }
  )
  out <- do.call(rbind, data)
  return(out)
}
