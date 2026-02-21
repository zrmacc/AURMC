# TabRMC() tests

test_that("TabRMC returns expected structure.", {
  df <- data.frame(
    idx = c(1, 1, 2, 2, 2),
    time = c(0, 1, 0, 1, 2),
    status = c(1, 0, 1, 1, 0),
    value = c(1, 1, 1, 1, 1)
  )
  out <- TabRMC(df)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("time", "nar", "y", "haz", "surv", "d", "exp") %in% names(out)))
})

test_that("TabRMC matches EstimatorR curve.", {
  df <- data.frame(
    idx = c(1, 1, 2, 2, 2),
    time = c(0, 1, 0, 1, 2),
    status = c(1, 0, 1, 1, 0),
    value = c(1, 1, 1, 1, 1)
  )
  tab <- TabRMC(df)
  est <- AURMC:::EstimatorR(
    idx = df$idx,
    status = df$status,
    time = df$time,
    value = df$value,
    return_auc = FALSE
  )
  expect_equal(tab$time, est$time)
  expect_equal(tab$exp, est$exp)
})
