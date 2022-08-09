test_that("Test value tabulation.", {
  
  # Case 1.
  df <- data.frame(
    idx = c(1, 1),
    time = c(1, 2),
    value = c(1, NA)
  )
  
  # Observed.
  obs <- ValueMatrixR(
    eval_times = c(0, 1, 2, 3),
    idx = df$idx,
    time = df$time,
    value = df$value
  )
  
  # Expected.
  exp <- c(0, 1, 1, 0)
  expect_equal(as.numeric(obs), exp)
  
  # Case 2.
  df <- data.frame(
    idx = c(1, 1),
    time = c(1, 2),
    value = c(1, 2)
  )
  
  # Observed.
  obs <- ValueMatrixR(
    eval_times = c(0, 1, 2, 3),
    idx = df$idx,
    time = df$time,
    value = df$value
  )
  
  # Expected.
  exp <- c(0, 1, 2, 0)
  expect_equal(as.numeric(obs), exp)
  
  
  # Case 3.
  df <- data.frame(
    idx = c(1, 1, 2, 2, 2, 3, 3, 3, 3),
    time = c(0, 0.5, 0, 1, 1.5, 0, 1, 2, 2.5),
    value = c(1, 2, 1, 2, NA, 1, 2, 3, NA)
  )
  
  # Observed.
  obs <- ValueMatrixR(
    eval_times = sort(unique(df$time)),
    idx = df$idx,
    time = df$time,
    value = df$value
  )
  
  # Expected.
  exp <- rbind(
    c(1, 2, 0, 0, 0, 0),
    c(1, 1, 2, 2, 0, 0),
    c(1, 1, 2, 2, 3, 3)
  )
  
  expect_equal(obs, exp)
})
