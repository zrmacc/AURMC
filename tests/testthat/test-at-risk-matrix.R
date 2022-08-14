test_that("Test at risk matrix.", {
  
  # Case 1.
  df <- data.frame(
    idx = c(1, 1),
    time = c(1, 2)
  )
  
  # Observed.
  obs <- AtRiskMatrixR(
    eval_times = c(1, 2, 3),
    idx = df$idx,
    time = df$time  
  )
  
  # Expected.
  exp <- c(1, 1, 0)
  expect_equal(as.numeric(obs), exp)
  
  # Case 2.
  df <- data.frame(
    idx = c(1, 1, 2, 2, 3, 3),
    time = c(0, 2, 0, 0.5, 0, 3)
  )
  
  # Observed.
  obs <- AtRiskMatrixR(
    eval_times = c(0, 1, 2, 3),
    idx = df$idx,
    time = df$time
  )
  
  # Expected.
  exp <- rbind(
    c(1, 1, 1, 0),
    c(1, 0, 0, 0),
    c(1, 1, 1, 1)
  )
  expect_equal(obs, exp)

})
