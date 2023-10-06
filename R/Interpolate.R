# Purpose: Wrapper for interpolation function.
# Updated: 2022-08-20

#' Interpolate
#' 
#' Linearly interpolations between each subject's measurements.
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
  
  # Convert index to numeric.
  idx_replaced <- FALSE
  idx <- numeric_idx <- NULL
  if (is.factor(data$idx)) {
    data$numeric_idx <- as.numeric(data$idx)
    idx_map <- data %>% dplyr::select(idx, numeric_idx) %>% unique()
    data$idx <- data$numeric_idx
    data$numeric_idx <- NULL
    idx_replaced <- TRUE
  }
  
  # Union grid with each subject's last time.
  time <- last_time <- NULL
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
  if (idx_replaced) {
    interpolated$numeric_idx <- interpolated$idx
    interpolated$idx <- NULL
    interpolated <- interpolated %>%
      dplyr::inner_join(idx_map, by = "numeric_idx")
    interpolated$numeric_idx <- NULL
  }
  
  return(interpolated)
}
