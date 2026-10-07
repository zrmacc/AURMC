# Purpose: Main estimation and inference function.
# Updated: 2024-11-09


#' Difference
#'
#' @param one_sample Data.frame of one sample results.
#' @param alpha Type I error.
#' @return Data.frame for the difference.
#' @noRd
Diff <- function(one_sample, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  arm <- auc <- est <- se <- NULL
  out <- one_sample %>%
    dplyr::summarise(
      stat = "A1-A0",
      est = auc[arm == 1] - auc[arm == 0],
      se = sqrt(se[arm == 1]^2 + se[arm == 0]^2)
    ) %>%
    dplyr::mutate(
      lower = est - z * se,
      upper = est + z * se,
      p = if (se == 0) {
        if (est == 0) 1 else 0
      } else {
        2 * stats::pnorm(q = abs(est) / se, lower.tail = FALSE)
      }
    )
  return(out)
}


#' Ratio
#'
#' @param one_sample Data.frame of one sample results.
#' @param alpha Type I error.
#' @return Data.frame for the difference.
#' @noRd
Ratio <- function(one_sample, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  arm <- auc <- est <- log_se <- se <- NULL
  auc0 <- one_sample$auc[one_sample$arm == 0]
  auc1 <- one_sample$auc[one_sample$arm == 1]
  if (length(auc0) != 1 || length(auc1) != 1) {
    stop("Ratio inference requires one result from each treatment arm.", call. = FALSE)
  }
  if (!is.finite(auc0) || !is.finite(auc1) || auc0 <= 0 || auc1 <= 0) {
    return(data.frame(
      stat = "A1/A0",
      est = auc1 / auc0,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      p = NA_real_
    ))
  }
  out <- one_sample %>%
    dplyr::summarise(
      stat = "A1/A0",
      est = auc[arm == 1] / auc[arm == 0],
      log_se = sqrt(se[arm == 1]^2 / auc[arm == 1]^2 + se[arm == 0]^2 / auc[arm == 0]^2)
    ) %>% 
    dplyr::mutate(
      lower = est * exp(-z * log_se),
      upper = est * exp(+z * log_se),
      se = est * log_se,
      p = if (log_se == 0) {
        if (est == 1) 1 else 0
      } else {
        2 * stats::pnorm(q = abs(log(est)) / log_se, lower.tail = FALSE)
      }
    ) %>% 
    dplyr::select(-log_se)
  return(out)
}


#' Difference and Ratio
#' 
#' @param arm0 Data.frame of arm 0 results.
#' @param arm1 Data.frame of arm 1 results.
#' @param alpha Type I error.
#' @return Data.frame of contrasts.
#' @noRd
DiffRatio <- function(arm0, arm1, alpha = 0.05) {
  
  one_sample <- rbind(arm0, arm1)
  split_data <- split(one_sample, one_sample$method)
  out <- lapply(split_data, function(df) {
    method <- unique(df$method)
    diff <- Diff(df, alpha = alpha)
    ratio <- Ratio(df, alpha = alpha)
    diff_ratio <- rbind(diff, ratio) %>%
      tibble::add_column(.before = "stat", method = method)
  })
  names(out) <- NULL
  out <- do.call(rbind, out)
  return(out)
}


#' Compare Areas Under Repeated Measures Curves
#' 
#' @param data Data.frame.
#' @param alpha Type I error.
#' @param arm_name Name of the column containing treatment arm. Coded as 0 for
#'   reference, 1 for treatment.
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
#' @return An object of class `AURMC` containing arm-specific estimates
#'   and difference and ratio contrasts. Ratio inference is returned as
#'   `NA` when either arm-specific area is non-positive.
#' @examples
#' example_data <- data.frame(
#'   idx = rep(1:4, each = 2),
#'   arm = rep(c(0, 0, 1, 1), each = 2),
#'   time = rep(c(0, 1), 4),
#'   status = rep(c(1, 0), 4),
#'   value = rep(c(1, 1, 2, 2), each = 2)
#' )
#' CompareAURMCs(example_data, tau = 1)
#' @export
CompareAURMCs <- function(
  data,
  alpha = 0.05,
  arm_name = "arm",
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
  
  # Format input data.
  arm <- idx <- status <- time <- value <- NULL
  data <- data %>%
    dplyr::rename(
      idx = {{idx_name}},
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}},
      value = {{value_name}}
    )

  ValidateCoreInput(data, check_arm = TRUE)
  
  # Censor after last.
  if (censor_after_last) {
    data <- CensorAfterLast(data)
  }

  # Encode subject identifiers and order records chronologically.
  data <- PrepareEstimatorInput(data)
  
  # Truncation time.
  arm_max_time <- tapply(data$time, data$arm, max)
  if (is.null(tau)) {
    tau <- min(arm_max_time)
  }
  if (length(tau) != 1 || !is.finite(tau) || tau <= 0 || tau > min(arm_max_time)) {
    stop(
      "`tau` must be positive and no greater than the maximum follow-up time in either arm.",
      call. = FALSE
    )
  }
  
  # Check input.
  InputCheck(data, check_arm = TRUE)
  
  # Perform 1 sample inference.
  arm0 <- AURMC(
    data = data %>% dplyr::filter(arm == 0),
    alpha = alpha,
    censor_after_last = FALSE,
    int_method = int_method,
    perturbations = perturbations,
    random_state = random_state,
    tau = tau
  )
  arm0 <- arm0 %>%
    tibble::add_column(.after = "method", arm = 0)
  
  arm1 <- AURMC(
    data = data %>% dplyr::filter(arm == 1),
    alpha = alpha,
    censor_after_last = FALSE,
    int_method = int_method,
    perturbations = perturbations,
    random_state = random_state,
    tau = tau
  )
  arm1 <- arm1 %>%
    tibble::add_column(.after = "method", arm = 1)
  
  # Contrast 1 sample results.
  contrast <- DiffRatio(arm0, arm1, alpha = alpha)
  
  # Output.
  out <- methods::new(
    Class = "AURMC",
    Arm0 = arm0,
    Arm1 = arm1,
    Contrast = contrast
  )
  return(out)
}
