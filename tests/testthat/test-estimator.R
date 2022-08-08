test_that("Test estimator.", {
  
  df <- data.frame(
    idx = c(1, 1, 2, 2, 2, 3, 3, 3, 3),
    status = c(1, 0, 1, 1, 0, 1, 1, 1, 0),
    time = c(0, 0.5, 0, 1, 1.5, 0, 1, 2, 2.5),
    value = c(1, NA, 1, 1, NA, 1, 1, 1, NA)
  )
  
  # Observed.
  obs <- EstimatorR(
    idx = df$idx,
    status = df$status,
    time = df$time,
    value = df$value,
    return_auc = FALSE
  )
  
  # Expected.
  n_times <- length(unique(df$time))
  expect_equal(obs$exp_value, rep(1, n_times))
  
  # Observed.
  obs <- EstimatorR(
    idx = df$idx,
    status = df$status,
    time = df$time,
    value = df$value,
    return_auc = TRUE,
    trunc_time = 1.25
  )
  
  # Expected.
  expect_equal(obs, 1.25)
  
})
