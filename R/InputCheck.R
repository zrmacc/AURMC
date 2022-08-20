# Purpose: Check input.
# Updated: 2022-08-14


#' Check Subject
#' 
#' @param idx Subject index.
#' @param status Event status.
#' @param time Observation time.
#' @return Logical.
#' @noRd
CheckSubj <- function(idx, status, time) {
  
  idx <- unique(idx)
  has_time_zero <- (time[1] == 0)
  time_zero_status <- (status[1] == 1)
  failed <- FALSE
  
  if (!has_time_zero | !time_zero_status) {
    failed <- TRUE
    warning(paste0("Subject ", idx, " lacks a record with time = 0 and status = 1."))
  }
  
  obs_end <- (status == 0 | status == 2)
  any_obs_end <- any(obs_end)
  
  if(!any_obs_end) {
    failed <- TRUE
    warning(paste0("Subject ", idx, " has no observation terminating event (status = 0 or status = 2)."))
  }
  
  sum_obs_end <- sum(obs_end)
  if(sum_obs_end > 1) {
    failed <- TRUE
    warning(paste0("Subject ", idx, " has multiple observation terminating events (status = 0 or status = 2)."))
  }
  return(failed)  
}


#' Input Check
#' 
#' Check for proper input formatting.
#' 
#' @param data Data.frame.
#' @return None.
InputCheck <- function(data) {
  
  idx <- status <- time <- NULL
  check <- data %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      failed = CheckSubj(idx, status, time)
    )
  failed <- any(check$failed)
  
  # Check treatment arm, if present.
  if ("arm" %in% colnames(data)) {
    # Check coding.
    arm_levels <- sort(unique(data$arm))
    is_proper_coding <- all(arm_levels == c(0, 1))
    if (!is_proper_coding) {
      failed <- TRUE
      warning("Treatment arm is improperly coded. Expecting two levels, c(0, 1).")
    }
    
    # Check for non-unique index.
    idx0 <- unique(data$idx[data$arm == 0])
    idx1 <- unique(data$idx[data$arm == 1])
    idx_overlap <- intersect(idx0, idx1)
    if (length(idx_overlap) > 0) {
      failed <- TRUE
      msg <- "The following subject indices appear in both treatment arms:\n "
      msg <- paste0(msg, paste(idx_overlap, collapse = ", "))
      warning(msg)
    }
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
  
  split_data <- split(x = data, f = data$idx)
  formatted_data <- lapply(split_data, function(df) {
    
    obs_end <- (df$status == 0 | df$status == 2)
    any_obs_end <- any(obs_end)
    
    if (any_obs_end) {
      return(df)
    }
    
    # Add censoring record.
    last_row <- df[nrow(df), ]
    last_row$status <- 0
    last_row$time <- last_row$time + 1e-4
    df <- rbind(df, last_row)
    return(df)
    
  })
  
  out <- do.call(rbind, formatted_data)
  rownames(out) <- NULL
  return(out)
}
