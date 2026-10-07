# Purpose: Check input.
# Updated: 2022-08-20


#' Check Subject
#' 
#' @param idx Subject index.
#' @param status Event status.
#' @param time Observation time.
#' @return Logical.
#' @noRd
CheckSubj <- function(idx, status, time, require_end = TRUE) {
  
  idx <- unique(idx)
  has_baseline <- any(time == 0 & status == 1)
  failed <- FALSE

  if (!has_baseline) {
    failed <- TRUE
    warning(paste0("Subject ", idx, " lacks a record with time = 0 and status = 1."))
  }
  
  obs_end <- (status == 0 | status == 2)
  any_obs_end <- any(obs_end)
  
  if(require_end && !any_obs_end) {
    failed <- TRUE
    warning(paste0("Subject ", idx, " has no observation terminating event (status = 0 or status = 2)."))
  }
  
  sum_obs_end <- sum(obs_end)
  if(sum_obs_end > 1) {
    failed <- TRUE
    warning(paste0("Subject ", idx, " has multiple observation terminating events (status = 0 or status = 2)."))
  }

  if (sum_obs_end == 1 && time[obs_end] != max(time)) {
    failed <- TRUE
    warning(paste0(
      "Subject ", idx,
      " has records after the observation terminating event."
    ))
  }
  return(failed)  
}


#' Arm Check
#' 
#' Check formatting of treatment arm.
#' 
#' @param arm Treatment arm.
#' @param idx Subject index.
#' @return Logical.
#' @noRd
CheckArm <- function(arm, idx) {
  failed <- FALSE
  
  # Check coding.
  arm_levels <- sort(unique(arm))
  is_proper_coding <- length(arm_levels) == 2 && all(arm_levels == c(0, 1))
  if (!is_proper_coding) {
    failed <- TRUE
    warning("Treatment arm is improperly coded. Expecting two levels, c(0, 1).")
  }
  
  # Check for non-unique index.
  idx0 <- unique(idx[arm == 0])
  idx1 <- unique(idx[arm == 1])
  idx_overlap <- intersect(idx0, idx1)
  if (length(idx_overlap) > 0) {
    failed <- TRUE
    msg <- "The following subject indices appear in both treatment arms:\n "
    msg <- paste0(msg, paste(idx_overlap, collapse = ", "))
    warning(msg)
  }
  
  return(failed)
}


#' Validate Core Input
#'
#' Validate fields used by the compiled estimators before coercion.
#'
#' @param data Data.frame with standardized column names.
#' @param check_arm Check the treatment-arm field?
#' @return None.
#' @noRd
ValidateCoreInput <- function(data, check_arm = FALSE) {
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("`data` must be a non-empty data.frame.", call. = FALSE)
  }

  required <- c("idx", "status", "time", "value")
  if (check_arm) {
    required <- c(required, "arm")
  }
  missing_names <- setdiff(required, names(data))
  if (length(missing_names) > 0) {
    stop(
      "Missing required column(s): ", paste(missing_names, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (anyNA(data$idx)) {
    stop("Subject indices must not be missing.", call. = FALSE)
  }
  if (!is.numeric(data$time) || anyNA(data$time) || any(!is.finite(data$time))) {
    stop("Observation times must be finite, non-missing numeric values.", call. = FALSE)
  }
  if (any(data$time < 0)) {
    stop("Observation times must be non-negative.", call. = FALSE)
  }
  if (!is.numeric(data$status) || anyNA(data$status) ||
      any(!is.finite(data$status)) || any(!data$status %in% c(0, 1, 2))) {
    stop("Status must be numeric and coded as 0, 1, or 2.", call. = FALSE)
  }
  if (!is.numeric(data$value)) {
    stop("Measurement values must be numeric (missing values are allowed).", call. = FALSE)
  }
  if (any(!is.na(data$value) & !is.finite(data$value))) {
    stop("Measurement values must be finite or missing.", call. = FALSE)
  }

  if (check_arm &&
      (!is.numeric(data$arm) || anyNA(data$arm) ||
       any(!is.finite(data$arm)) || any(!data$arm %in% c(0, 1)))) {
    stop("Treatment arm must be numeric and coded as 0 or 1.", call. = FALSE)
  }

  invisible(NULL)
}


#' Prepare Estimator Input
#'
#' Encode arbitrary subject identifiers as consecutive integers and order each
#' subject's records chronologically, as required by the compiled routines.
#'
#' @param data Data.frame with standardized column names.
#' @return Prepared data.frame.
#' @noRd
PrepareEstimatorInput <- function(data) {
  data$idx <- match(data$idx, unique(data$idx))
  data <- data[order(data$idx, data$time), , drop = FALSE]
  rownames(data) <- NULL
  data
}


#' Input Check
#' 
#' Check for proper input formatting.
#' 
#' @param data Data.frame.
#' @param check_arm Check arm?
#' @return None.
#' @noRd
InputCheck <- function(data, check_arm = FALSE, require_end = TRUE) {
  ValidateCoreInput(data, check_arm = check_arm)
  
  idx <- status <- time <- NULL
  check <- data %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      failed = CheckSubj(idx, status, time, require_end = require_end),
      .groups = "drop"
    )
  failed <- any(check$failed)
  
  # Check treatment arm.
  if (check_arm) {
    failed <- any(failed, CheckArm(data$arm, data$idx))
  }
  
  if (failed) {
    stop("Input check failed.")
  }

  return(invisible(NULL))
}


#' Censor After Last
#' 
#' Introduce a censoring after the last event if no observation
#' terminating event is present.
#' 
#' @param data Data.frame.
#' @return None.
CensorAfterLast <- function(data) {
  ValidateCoreInput(data, check_arm = "arm" %in% names(data))
  data <- data[order(data$idx, data$time), , drop = FALSE]
  rownames(data) <- NULL
  
  split_data <- split(x = data, f = data$idx)
  formatted_data <- lapply(split_data, function(df) {
    
    obs_end <- (df$status == 0 | df$status == 2)
    any_obs_end <- any(obs_end)
    
    if (any_obs_end) {
      return(df)
    }
    
    # Add censoring record.
    last_row <- df[which.max(df$time), , drop = FALSE]
    last_row$status <- 0
    last_row$time <- last_row$time + 1e-4
    df <- rbind(df, last_row)
    return(df)
    
  })
  
  out <- do.call(rbind, formatted_data)
  rownames(out) <- NULL
  return(out)
}
