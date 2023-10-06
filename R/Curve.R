# Purpose: Plot the integrand of the AUC.
# Updated: 2022-08-20


#' Estimate Curve
#' 
#' @param data Data.frame.
#' @param censor_after_last Introduce censoring after the last event *if* no
#'   observation-terminating event is present.
#' @param idx_name Name of column containing a unique subject index.
#' @param status_name Name of column containing the status. Must be coded as 0
#'   for censoring, 1 for a measurement, 2 for death. Each subject should have
#'   an observation-terminating event, either censoring or death.
#' @param time_name Name of column containing the observation time.
#' @param value_name Name of the column containing the measurement.
#' @return Function.
#' @importFrom dplyr "%>%"
#' @export
Curve <- function(
  data,
  censor_after_last = TRUE,
  idx_name = "idx",
  status_name = "status",
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
  
  # Convert index to numeric.
  if (is.factor(data$idx)) {
    data$idx <- as.numeric(data$idx)
  }
  
  # Censor after last.
  if (censor_after_last) {
    data <- CensorAfterLast(data)
  }
  
  # Check input.
  InputCheck(data, check_arm = FALSE)
  
  # Tabulate repeated measures curve. 
  est <- EstimatorR(
    idx = data$idx,
    status = data$status,
    time = data$time,
    value = data$value,
    return_auc = FALSE
  )
  
  # Integrand.
  fn <- stats::approxfun(
    x = est$time,
    y = est$exp,
    rule = 2,
    f = 0.5
  )
  return(fn)
}