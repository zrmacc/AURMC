test_that("Test interpolation.", {
  
  data <- data.frame(
    idx = c(1, 1, 2, 2),
    status = c(1, 0, 1, 2),
    time = c(0, 1, 0, 1),
    value = c(0, 1, 0, -1)
  )
  
  obs <- InterpolateR(
    idx = data$idx,
    status = data$status,
    time = data$time,
    value = data$value,
    n_points = 4
  )
  
  seq_0_1_4 <- seq(from = 0, to = 1, length.out = 4)
  exp <- data.frame(
    idx = c(rep(1, 4), rep(2, 4)),
    status = c(1, 1, 1, 0, 1, 1, 1, 2),
    time = rep(seq_0_1_4, 2),
    value = c(seq_0_1_4, -1 * seq_0_1_4)
  )
  expect_equal(obs, exp)
  
})