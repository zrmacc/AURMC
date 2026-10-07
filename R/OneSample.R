# Purpose: Main estimation and inference function.
# Updated: 2024-11-09


#' Area Under the Repeated Measures Curve
#' 
#' @param data Data.frame.
#' @param alpha Type I error.
#' @param censor_after_last Introduce censoring after the last event *if* no
#'   observation-terminating event is present.
#' @param idx_name Name of column containing a unique subject index.
#' @param int_method Integration method, selected from "left", "right", "trapezoid".
#' @param perturbations Number of perturbations to use for bootstrap inference.
#'   If \code{NULL}, only analytical inference is performed.
#' @param random_state Seed to ensure perturbations are reproducible. 
#' @param status_name Name of column containing the status. Must be coded as 0
#'   for censoring, 1 for a measurement, 2 for death. Each subject should have
#'   an observation-terminating event, either censoring or death.
#' @param tau Truncation time.
#' @param time_name Name of column containing the observation time.
#' @param value_name Name of the column containing the measurement.
#' @return A data.frame containing the method, truncation time, area estimate,
#'   standard error, confidence limits, and p-value.
#' @examples
#' example_data <- data.frame(
#'   idx = rep(1:3, each = 2),
#'   time = rep(c(0, 1), 3),
#'   status = rep(c(1, 0), 3),
#'   value = c(1, 1, 2, 2, 3, 3)
#' )
#' AURMC(example_data, tau = 1)
#' @export
AURMC <- function(
  data,
  alpha = 0.05,
  censor_after_last = TRUE,
  idx_name = "idx",
  int_method = "trapezoid",
  perturbations = NULL,
  random_state = 0,
  status_name = "status",
  tau = NULL,
  time_name = "time",
  value_name = "value"
) {
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  int_method <- match.arg(int_method, c("left", "right", "trapezoid"))
  if (!is.null(perturbations) &&
      (length(perturbations) != 1 || !is.finite(perturbations) ||
       perturbations < 2 || perturbations != as.integer(perturbations))) {
    stop("`perturbations` must be NULL or an integer of at least 2.", call. = FALSE)
  }
  
  # Format input data.
  data <- data %>%
    dplyr::rename(
      idx = {{idx_name}},
      status = {{status_name}},
      time = {{time_name}},
      value = {{value_name}}
    )

  ValidateCoreInput(data, check_arm = FALSE)
  
  # Censor after last.
  if (censor_after_last) {
    data <- CensorAfterLast(data)
  }

  # Encode subject identifiers and order records chronologically.
  data <- PrepareEstimatorInput(data)
  
  # Truncation time.
  if (is.null(tau)) {
    tau <- max(data$time)
  }
  if (length(tau) != 1 || !is.finite(tau) || tau <= 0 || tau > max(data$time)) {
    stop(
      "`tau` must be a single positive number no greater than the maximum follow-up time.",
      call. = FALSE
    )
  }
  
  # Check input.
  InputCheck(data, check_arm = FALSE)
  
  # Calculate AUC.
  auc <- EstimatorR(
    idx = data$idx,
    status = data$status,
    time = data$time,
    value = data$value,
    int_method = int_method,
    trunc_time = tau,
    return_auc = TRUE
  )
  
  # Calculate SE.
  n <- length(unique(data$idx))
  psi <- InfluenceR(
    idx = data$idx,
    int_method = int_method,
    status = data$status,
    time = data$time,
    trunc_time = tau,
    value = data$value
  )
  se <- sqrt(mean(psi$psi^2) / n)
  
  # Asymptotic output.
  z <- stats::qnorm(p = 1 - alpha / 2)
  if (isTRUE(all.equal(se, 0))) {
    p <- if (isTRUE(all.equal(auc, 0))) 1 else 0
  } else {
    p <- stats::pchisq(q = (auc / se)^2, df = 1, lower.tail = FALSE)
  }
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
      int_method = int_method,
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
    lower_tail <- (sum(auc_jitter <= 0) + 1) / (perturbations + 1)
    upper_tail <- (sum(auc_jitter >= 0) + 1) / (perturbations + 1)
    boot_p <- min(1, 2 * min(lower_tail, upper_tail))
    
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


#' Tabulate the Repeated Measures Curve
#' 
#' @param data Data.frame.
#' @param censor_after_last Introduce censoring after the last event *if* no
#'   observation-terminating event is present.
#' @param idx_name Name of column containing a unique subject index.
#' @param status_name Name of column containing the status. Must be coded as 0
#'   for censoring, 1 for a measurement, 2 for death. Each subject should have
#'   an observation-terminating event, either censoring or death.
#' @param tau Truncation time.
#' @param time_name Name of column containing the observation time.
#' @param value_name Name of the column containing the measurement.
#' @return A data.frame tabulating time, number and proportion at risk,
#'   terminal-event hazard, left-limit survival \eqn{\hat S(t-)}, observed-value
#'   mean, and the estimated repeated measures curve.
#' @export
TabRMC <- function(
    data,
    censor_after_last = TRUE,
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

  ValidateCoreInput(data, check_arm = FALSE)
  
  # Censor after last.
  if (censor_after_last) {
    data <- CensorAfterLast(data)
  }

  # Encode subject identifiers and order records chronologically.
  data <- PrepareEstimatorInput(data)
  
  # Truncation time.
  if (is.null(tau)) {
    tau <- max(data$time)
  }
  if (length(tau) != 1 || !is.finite(tau) || tau <= 0 || tau > max(data$time)) {
    stop(
      "`tau` must be a single positive number no greater than the maximum follow-up time.",
      call. = FALSE
    )
  }
  
  # Check input.
  InputCheck(data, check_arm = FALSE)
  
  # Calculate AUC.
  out <- EstimatorR(
    idx = data$idx,
    status = data$status,
    time = data$time,
    value = data$value,
    trunc_time = tau,
    return_auc = FALSE
  )
  return(out)
}
