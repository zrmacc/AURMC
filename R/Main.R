# Purpose: Main estimation and inference function.
# Updated: 2022-08-14

#' Area Under the Repeated Measures Curve
#' 
#' @param data Data.frame.
#' @param alpha Type I error.
#' @param idx_name Name of column containing a unique subject index.
#' @param status_name Name of column containing the status. Must be coded as 0 for
#'   censoring, 1 for a measurement, 2 for death. Each subject should have an 
#'   observation-terminating event, either censoring or death.
#' @param tau Truncation time.
#' @param time_name Name of column containing the observation time.
#' @param value_name Name of the column containing the measurement.
#' @return Data.frame.
#' @importFrom dplyr "%>%"
#' @export 
AURMC <- function(
  data,
  alpha = 0.05,
  idx_name = "idx",
  status_name = "status",
  tau = NULL,
  time_name = "time",
  value_name = "value"
) {
  
  # Format input data.
  data <- data %>%
    dplyr::rename(
      idx = {{idx_name}},
      status = {{status_name}},
      time = {{time_name}},
      value = {{value_name}}
    )
  
  # Truncation time.
  if (is.null(tau)) {
    tau <- max(data$time)
  }
  
  # Check input.
  InputCheck(data)
  
  # Calculate AUC.
  auc <- EstimatorR(
    idx = data$idx,
    status = data$status,
    time = data$time,
    value = data$value,
    trunc_time = tau,
    return_auc = TRUE
  )
  
  # Calculate SE.
  n <- length(unique(data$idx))
  psi <- InfluenceR(
    idx = data$idx,
    status = data$status,
    time = data$time,
    trunc_time = tau,
    value = data$value
  )
  se <- sqrt(mean(psi$psi^2) / n)
  
  # Prepare output.
  z <- stats::qnorm(p = 1 - alpha / 2)
  p <- stats::pchisq(q = (auc / se)^2, df = 1, lower.tail = FALSE)
  out <- data.frame(
    tau = tau,
    auc = auc,
    se = se
  )
  out$lower <- out$auc - z * out$se
  out$upper <- out$auc + z * out$se
  out$p <- p
  return(out)
}