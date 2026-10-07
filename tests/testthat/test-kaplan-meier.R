test_that("Kaplan-Meier tabulation returns left-limit survival.", {

  df <- data.frame(
    idx = c(1, 2, 3, 4),
    time = c(1, 2, 3, 4),
    status = c(0, 2, 0, 2)
  )
  
  obs <- KaplanMeierR(
    eval_times = seq(from = 0, to = 5),
    idx = df$idx,
    status = df$status,
    time = df$time
  )
  
  expect_equal(obs$nar, c(4, 4, 3, 2, 1, 0))
  expect_equal(obs$haz, c(0, 0, 1/3, 0, 1, NaN))
  expect_equal(obs$surv, c(1, 1, 1, 2/3, 2/3, 0))

  # At a death time, surv is S(t-): the current hazard jump is not included.
  # This discriminates the required convention from post-jump S(t).
  expect_equal(obs$surv[obs$time == 2], 1)
  expect_equal(obs$surv[obs$time == 4], 2/3)

})
