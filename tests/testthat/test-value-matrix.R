test_that("Test value tabulation.", {
  
  df <- data.frame(
    idx = c(1, 1, 2, 2, 2, 3, 3, 3, 3),
    status = c(1, 2, 1, 1, 0, 1, 1, 1, 2),
    time = c(0, 0.5, 0, 1, 1.5, 0, 1, 2, 2.5),
    value = c(1, NA, 1, 2, NA, 1, 2, 3, NA)
  )
  
  # Observed.
  obs <- ValueMatrixR(
    eval_times = sort(unique(df$time)),
    idx = df$idx,
    status = df$status,
    time = df$time,
    value = df$value
  )
  
  # Expected.
  exp <- rbind(
    c(1, 0, 0, 0, 0, 0),
    c(1, 1, 2, 0, 0, 0),
    c(1, 1, 2, 2, 3, 0)
  )
  
  expect_equal(obs, exp)
})
