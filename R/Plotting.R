# Purpose: Plotting.
# Updated: 2023-10-30


#' Plot Two Sample
#' 
#' @param data Data.frame.
#' @param arm_name Name of the column containing treatment arm. Coded as 0 for
#'   reference, 1 for treatment.
#' @param censor_after_last Introduce censoring after the last event *if* no
#'   observation-terminating event is present.
#' @param color_labs Color labels.
#' @param colors Colors for control and treatment respectively.
#' @param eval_pts Evaluation points between 0 and tau.
#' @param idx_name Name of column containing a unique subject index.
#' @param status_name Name of column containing the status. Must be coded as 0
#'   for censoring, 1 for a measurement, 2 for death. Each subject should have
#'   an observation-terminating event, either censoring or death.
#' @param show_auc Show area under the curve?
#' @param show_nar Show number at risk?
#' @param tau Truncation time.
#' @param time_name Name of column containing the observation time.
#' @param title Plot title.
#' @param value_name Name of the column containing the measurement.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis labels.
#' @param x_name X-axis name.
#' @param x_max X-axis upper limit.
#' @param y_name Y-axis name.
#' @param y_lim Y-axis limits.
#' @return ggplot
#' @export
PlotTwoSample <- function(
  data,
  arm_name = "arm",
  censor_after_last = TRUE,
  color_labs = c("Ctrl", "Trt"),
  colors = c("#EFC000FF", "#6385B8"),
  eval_pts = 200,
  idx_name = "idx",
  show_auc = FALSE,
  show_nar = TRUE,
  status_name = "status",
  tau = NULL,
  time_name = "time",
  title = NULL,
  value_name = "value",
  x_breaks = NULL,
  x_labs = NULL,
  x_name = "Time",
  x_max = NULL,
  y_name = "Expected Value",
  y_lim = NULL
) {
  
  # Format input data.
  arm <- idx <- status <- time <- value <- NULL
  data <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
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
  
  # Truncation time.
  if (is.null(x_max)) {
    x_max <- max(data$time)
  }
  
  if (is.null(tau)) {
    tau <- x_max
  }
  
  # Time axis.
  if (is.null(x_breaks)) {
    x_breaks <- seq(from = 0, to = tau, by = tau / 10)
  }
  
  if (is.null(x_labs)) {
    x_labs <- x_breaks
  }
  
  # Plotting frame.
  eval_times <- seq(from = 0, to = tau, length.out = eval_pts)
  eval_times <- sort(union(eval_times, x_breaks))
  
  df <- lapply(c(0, 1), function(a) {
    df <- data %>%
      dplyr::filter(arm == a)
    df <- EstimatorR(
      idx = df$idx,
      status = df$status,
      time = df$time,
      value = df$value,
      return_auc = FALSE,
      eval_times = eval_times,
      trunc_time = tau
    )  
    df$arm <- a
    return(df)
  })
  df <- do.call(rbind, df)
  df$arm <- factor(
    df$arm,
    levels = c(0, 1),
    labels = color_labs
  )

  # Construct plot.
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    ) + 
    ggplot2::geom_step(
      data = df, 
      ggplot2::aes(x = time, y = exp, color = arm),
      linewidth = 1
    )
  
  if (show_auc) {
    q <- q + 
      ggplot2::geom_ribbon(
        data = df,
        ggplot2::aes(x = time, ymin = 0, ymax = exp, fill = arm),
        alpha = 0.2
      )
  }
  
  # Options.
  q <- q + 
    ggplot2::scale_color_manual(
      name = NULL, 
      labels = color_labs, 
      values = colors
    ) + ggplot2::scale_fill_manual(
      name = NULL, 
      labels = color_labs, 
      values = colors
    ) +
    ggplot2::scale_x_continuous(
      name = x_name,
      breaks = x_breaks,
      labels = x_labs,
      limits = c(0, x_max)
    ) +
    ggplot2::scale_y_continuous(
      name = y_name,
      limits = y_lim
    ) + 
    ggplot2::ggtitle(
      label = title
    )
  
  if (show_nar) {
    
    df_nar <- df %>%
      dplyr::filter(time %in% x_breaks)
    
    nar <- NULL
    q_nar <- ggplot2::ggplot(data = df_nar) +
      ggplot2::theme_bw() + 
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = time, y = arm, label = nar)
      ) +
      ggplot2::scale_x_continuous(
        breaks = x_breaks,
        name = x_name,
        labels = x_labs,
        limits = c(0, x_max)
      ) + 
      ggplot2::scale_y_discrete(
        name = NULL
      )
    
    q_final <- cowplot::plot_grid(
      plotlist = list(q, q_nar),
      align = "v",
      axis = "l",
      ncol = 1,
      rel_heights = c(3, 1)
    )
  } else {
    
    q_final <- q
    
  }
  return(q_final)
}

