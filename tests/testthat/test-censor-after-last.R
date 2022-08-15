test_that("Test censor after last.", {
  
  # Case 1.
  df <- data.frame(
    idx = c(1, 1),
    status = c(1, 0),
    time = c(0, 1),
    value = c(1, 1)
  )
  
  # Observed.
  obs <- CensorAfterLast(df)
  expect_equal(obs, df)
  
  # Case 2.
  df <- data.frame(
    idx = c(1, 1),
    status = c(1, 2),
    time = c(0, 1),
    value = c(1, 1)
  )
  
  # Observed.
  obs <- CensorAfterLast(df)
  expect_equal(obs, df)
  
  # Case 3.
  df <- data.frame(
    idx = c(1, 1),
    status = c(1, 1),
    time = c(0, 1),
    value = c(1, 1)
  )
  
  # Observed.
  obs <- CensorAfterLast(df)
  exp <- rbind(df, c(1, 0, 1 + 1e-4, 1))
  expect_equal(obs, exp)

})
