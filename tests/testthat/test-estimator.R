test_that("Test estimator.", {
  
  df <- data.frame(
    idx = c(1, 1, 2, 2, 2, 3, 3, 3, 3),
    status = c(1, 0, 1, 1, 0, 1, 1, 1, 0),
    time = c(0, 0.5, 0, 1, 1.5, 0, 1, 2, 2.5),
    value = c(1, NA, 1, 1, NA, 1, 1, 1, NA)
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
  n_times <- length(unique(df$time))
  expect_equal(obs$exp, rep(1, n_times))
  
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


test_that("Test integration method.", {
  
  df <- data.frame(
    idx = c(1, 1, 1, 1, 1, 1),
    time = c(0, 1, 2, 3, 4, 5),
    status = c(1, 1, 1, 1, 1, 0),
    value = c(5, 4, 3, 2, 1, 0)
  )
  
  # Left sum.
  obs <- AURMC(df, int_method = "left")
  exp <- 5 + 4 + 3 + 2 + 1
  expect_equal(obs$auc[1], exp)
  
  # Right sum.
  obs <- AURMC(df, int_method = "right")
  exp <- 4 + 3 + 2 + 1 + 0
  expect_equal(obs$auc[1], exp)
  
  # Trapezoidal sum.
  obs <- AURMC(df, int_method = "trapezoid")
  exp <- 0.5 * (5 + 2 * 4 + 2 * 3 + 2 * 2 + 2 * 1 + 0)
  expect_equal(obs$auc[1], exp)
  
})
