# Purpose: Simulate data.
# Updated: 2022-08-06


#' Simulate Subject
#' 
#' @param censoring_rate Rate for the time to death.
#' @param death_rate Rate for the time to death.
#' @param idx Subject index.
#' @param tau Truncation time.
#' @return Data.frame.
#' @noRd
SimSubj <- function(
    censoring_rate,
    death_rate,
    idx,
    tau
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
  obs_values <- stats::rnorm(n = length(obs_times))
  
  # Output.
  out <- data.frame(
    idx = idx,
    time = c(obs_times, time),
    status = c(rep(1, length(obs_times)), status),
    values = c(obs_values, NA)
  )
  return(out)
}



#' Simulate Repeated Measures Data
#' 
#' @param censoring_rate Rate for the time to death.
#' @param death_rate Rate for the time to death.
#' @param n Number of subjects.
#' @param tau Truncation time.
#' @return Data.frame.
#' @export
GenData <- function(
    censoring_rate,
    death_rate,
    n,
    tau
) {
  
  data <- lapply(seq_len(n), function(i) {
    SimSubj(censoring_rate, death_rate, i, tau)
  })
  out <- do.call(rbind, data)
  return(out)
}


