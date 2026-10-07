test_that("Test interpolation.", {
  
  # Case 1.
  data <- data.frame(
    idx = c(1, 1, 2, 2),
    status = c(1, 0, 1, 2),
    time = c(0, 1, 0, 1),
    value = c(0, 1, 0, -1)
  )
  
  grid <- c(0, 0.5, 1.0)
  obs <- Interpolate(
    data = data,
    grid = grid
  )
  
  exp <- data.frame(
    idx = c(rep(1, 3), rep(2, 3)),
    status = c(1, 1, 0, 1, 1, 2),
    time = rep(grid, 2),
    value = c(grid, -1 * grid)
  )
  expect_equal(obs, exp)
  
  # Case 2.
  data <- data.frame(
    idx = c(1, 1),
    status = c(1, 2),
    time = c(0, 1),
    value = c(0, 1)
  )
  
  grid <- c(0, 1, 2)
  obs <- Interpolate(
    data = data,
    grid = grid,
    rm_na = FALSE
  )
  
  exp <- data.frame(
    idx = c(1, 1, 1),
    status = c(1, 2, 2),
    time = c(0, 1, 2),
    value = c(0, 1, NA)
  )
  expect_equal(obs, exp)
  
})

test_that("Interpolate preserves character identifiers and orders records.", {
  data <- data.frame(
    idx = c("b", "a", "b", "a"),
    status = c(2, 0, 1, 1),
    time = c(1, 1, 0, 0),
    value = c(-1, 1, 0, 0)
  )
  out <- Interpolate(data, grid = c(0, 0.5, 1))
  expect_type(out$idx, "character")
  expect_setequal(unique(out$idx), c("a", "b"))
  expect_equal(out$value[out$idx == "a"], c(0, 0.5, 1))
  expect_equal(out$value[out$idx == "b"], c(0, -0.5, -1))
})
