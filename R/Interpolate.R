# Purpose: Wrapper for interpolation function.
# Updated: 2022-08-20

#' Interpolate
#' 
#' Linearly interpolates between each subject's measurements.
#' The input data should contain no missing values. 
#'
#' @section Notes:
#' The grid of points will be augmented to include the time of each
#' subjects observation terminating event.
#'  
#' @param data Data.frame.
#' @param grid Grid of unique points at which to interpolate.
#' @param idx_name Name of column containing a unique subject index.
#' @param status_name Name of column containing the status. Must be coded as 0
#'   for censoring, 1 for a measurement, 2 for death. Each subject should have
#'   an observation-terminating event, either censoring or death.
#' @param rm_na Remove records interpolated to NA?
#' @param time_name Name of column containing the observation time.
#' @param value_name Name of the column containing the measurement.
#' @return A data.frame containing subject index, status, time, and the
#'   interpolated value.
#' @examples
#' example_data <- data.frame(
#'   idx = c(1, 1, 2, 2),
#'   status = c(1, 0, 1, 2),
#'   time = c(0, 1, 0, 1),
#'   value = c(0, 1, 0, -1)
#' )
#' Interpolate(example_data, grid = c(0, 0.5, 1))
#' @export 
Interpolate <- function(
  data,  
  grid, 
  idx_name = "idx",
  status_name = "status",
  rm_na = TRUE,
  time_name = "time",
  value_name = "value"
) {
  if (!is.numeric(grid) || length(grid) == 0 || anyNA(grid) ||
      any(!is.finite(grid)) || any(grid < 0)) {
    stop("`grid` must contain finite, non-negative numeric values.", call. = FALSE)
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
  
  # Encode subject identifiers and order records chronologically.
  original_idx <- unique(data$idx)
  data <- PrepareEstimatorInput(data)
  InputCheck(data, check_arm = FALSE)
  
  # Union grid with each subject's last time.
  idx <- time <- last_time <- NULL
  last_times <- data %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      last_time = max(time)
    ) %>%
    dplyr::pull(last_time)
  grid <- c(grid, last_times)
  grid <- sort(unique(grid))
  
  # Interpolate to grid.
  interpolated <- InterpolateR(
    grid = grid,
    idx = data$idx,
    status = data$status,
    time = data$time,
    value = data$value
  )
  
  # Remove NA.
  if (rm_na) {
    value <- NULL
    interpolated <- interpolated %>%
      dplyr::filter(!is.na(value))
  }
  
  # Restore original index.
  interpolated$idx <- original_idx[as.integer(interpolated$idx)]
  
  return(interpolated)
}
