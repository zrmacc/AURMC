test_that("Test curve tabulation.", {
  
  df <- data.frame(
    idx = c(1, 1, 1, 2, 2, 3, 3, 3),
    status = c(1, 1, 0, 1, 2, 1, 1, 0),
    time = c(0, 1, 2, 0, 1, 0, 1, 2),
    value = c(3, 2, 1, 1, 2, 2, 3, 1)
  )
  expect_error(InputCheck(df), NA)
  
  # Observed.
  obs <- EstimatorR(
    idx = df$idx,
    status = df$status,
    time = df$time,
    value = df$value,
    return_auc = FALSE
  )
  
  # Expected.
  df_km <- df %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      t_max = max(time),
      s_max = status[time == t_max]
    ) %>% dplyr::rename(time = t_max, status = s_max)
  
  km_tab <- KaplanMeierR(
    eval_times = c(0, 1, 2),
    idx = df_km$idx,
    status = df_km$status,
    time = df_km$time
  )
  
  df_sum <- df %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      d = sum(value),
      y = dplyr::n()
    )
  
  # Checking calculation of integrand. 
  exp <- km_tab$surv * (df_sum$d / df_sum$y)
  expect_equal(obs$exp, exp)
  
  # Checking curve.
  curve <- Curve(df)
  expect_equal(curve(obs$time), exp)
  
})
