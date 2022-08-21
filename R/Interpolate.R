# Purpose: Wrapper for interpolation function.
# Updated: 2022-08-20

#' Interpolate
#' 
#' Linearly interpolations between each subject's measurements.
#' The input data should contain no missing values. 
#'
#' @section Notes:
#' The grid of points will be agumented to include the time of each
#' subjects observating terminating event.
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
#' @return Data.frame.
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
  
  # Format input data.
  data <- data %>%
    dplyr::rename(
      idx = {{idx_name}},
      status = {{status_name}},
      time = {{time_name}},
      value = {{value_name}}
    )
  
  # Union grid with each subject's last time.
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
    interpolated <- interpolated %>%
      dplyr::filter(!is.na(value))
  }
  
  return(interpolated)
}