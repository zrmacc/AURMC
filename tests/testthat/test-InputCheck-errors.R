# InputCheck() error cases (warnings expected before stop)

test_that("InputCheck errors when subject lacks time 0 and status 1.", {
  df <- data.frame(
    idx = c(1, 1),
    time = c(0.5, 1),
    status = c(1, 0),
    value = c(1, 1)
  )
  expect_error(suppressWarnings(InputCheck(df, check_arm = FALSE)), "Input check failed")
})

test_that("InputCheck errors when subject has no terminating event.", {
  df <- data.frame(
    idx = c(1, 1),
    time = c(0, 1),
    status = c(1, 1),
    value = c(1, 1)
  )
  expect_error(suppressWarnings(InputCheck(df, check_arm = FALSE)), "Input check failed")
})

test_that("InputCheck errors when arm is not 0/1.", {
  df <- data.frame(
    idx = c(1, 1, 2, 2),
    arm = c(0, 0, 2, 2),
    time = c(0, 1, 0, 1),
    status = c(1, 0, 1, 0),
    value = c(1, 1, 1, 1)
  )
  expect_error(suppressWarnings(InputCheck(df, check_arm = TRUE)), "Input check failed")
})

test_that("InputCheck errors when same subject in both arms.", {
  df <- data.frame(
    idx = c(1, 1, 1, 1),
    arm = c(0, 0, 1, 1),
    time = c(0, 1, 0, 1),
    status = c(1, 0, 1, 0),
    value = c(1, 1, 1, 1)
  )
  expect_error(suppressWarnings(InputCheck(df, check_arm = TRUE)), "Input check failed")
})
