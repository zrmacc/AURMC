# GenData() tests

test_that("GenData returns expected columns.", {
  withr::local_seed(42)
  out <- GenData(censoring_rate = 0.5, death_rate = 0.25, n = 20, tau = 5)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("idx", "time", "status", "value"))
  expect_true(all(out$status %in% c(0, 1, 2)))
  expect_true(all(out$time >= 0))
  expect_true(all(out$idx >= 1))
})

test_that("GenData with last_missing produces NA at last value.", {
  withr::local_seed(42)
  out <- GenData(
    censoring_rate = 0.5, death_rate = 0.25, n = 10, tau = 5,
    last_missing = TRUE
  )
  # Last row per subject should have NA value when last_missing is TRUE
  last_per_subj <- out %>%
    dplyr::group_by(idx) %>%
    dplyr::slice_tail(n = 1)
  expect_true(any(is.na(last_per_subj$value)))
})

test_that("GenData passes InputCheck when censor_after_last applied.", {
  withr::local_seed(123)
  out <- GenData(censoring_rate = 0.3, death_rate = 0.2, n = 15, tau = 5)
  # After censoring after last, data should pass check (each subject has baseline and terminating event)
  out_censored <- CensorAfterLast(out)
  expect_error(InputCheck(out_censored, check_arm = FALSE), NA)
})

test_that("GenData validates simulation arguments.", {
  expect_error(GenData(-1, 0.25, 10, 2), "Event rates")
  expect_error(GenData(0.5, 0.25, 0, 2), "positive integer")
  expect_error(GenData(0.5, 0.25, 10, 0), "positive number")
  expect_error(GenData(0.5, 0.25, 10, 2, last_missing = NA), "TRUE or FALSE")
})
