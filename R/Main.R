# Purpose: Main estimation and inference function.
# Updated: 2022-08-20


#' Area Under the Repeated Measures Curve
#' 
#' @param data Data.frame.
#' @param alpha Type I error.
#' @param censor_after_last Introduce censoring after the last event *if* no
#'   observation-terminating event is present.
#' @param idx_name Name of column containing a unique subject index.
#' @param perturbations Number of perturbations to use for bootstrap inference.
#'   If \code{NULL}, only analytical inference is performed.
#' @param random_state Seed to ensure perturbations are reproducible. 
#' @param status_name Name of column containing the status. Must be coded as 0
#'   for censoring, 1 for a measurement, 2 for death. Each subject should have
#'   an observation-terminating event, either censoring or death.
#' @param tau Truncation time.
#' @param time_name Name of column containing the observation time.
#' @param value_name Name of the column containing the measurement.
#' @return Data.frame.
#' @importFrom dplyr "%>%"
#' @export
AURMC <- function(
  data,
  alpha = 0.05,
  censor_after_last = TRUE,
  idx_name = "idx",
  perturbations = NULL,
  random_state = 0,
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
  
  # Censor after last.
  if (censor_after_last) {
    data <- CensorAfterLast(data)
  }
  
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
  
  # Asymptotic output.
  z <- stats::qnorm(p = 1 - alpha / 2)
  p <- stats::pchisq(q = (auc / se)^2, df = 1, lower.tail = FALSE)
  out <- data.frame(
    method = "asymptotic",
    tau = tau,
    auc = auc,
    se = se
  )
  out$lower <- out$auc - z * out$se
  out$upper <- out$auc + z * out$se
  out$p <- p
  
  # Run perturbation.
  if (!is.null(perturbations)) {
    set.seed(random_state)
    deltas <- PerturbationR(
      idx = data$idx,
      perturbations = perturbations,
      status = data$status,
      time = data$time,
      trunc_time = tau,
      value = data$value
    )
    
    # Bootstrap SE.
    boot_se <- stats::sd(deltas)
    
    # Bootstrap CI.
    auc_jitter <- auc + deltas
    boot_ci <- stats::quantile(auc_jitter, probs = c(alpha / 2, 1 - alpha / 2))
    boot_ci <- as.numeric(boot_ci)
    
    # Bootstrap P.
    boot_p <- 2 * mean(sign(auc_jitter) != sign(auc))
    boot_p <- max(boot_p, 1 / perturbations)
    boot_p <- min(boot_p, 1)
    
    # Bootstrap results.
    out_boot <- data.frame(
      method = "bootstrap",
      tau = tau,
      auc = auc,
      se = boot_se,
      lower = boot_ci[1],
      upper = boot_ci[2],
      p = boot_p
    )
    out <- rbind(out, out_boot)
  }
  return(out)
}
