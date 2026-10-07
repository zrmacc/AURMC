# CompareAURMCs() tests

test_that("CompareAURMCs returns AURMC object with expected slots.", {
  df <- data.frame(
    idx = c(1, 1, 2, 2, 3, 3, 4, 4),
    arm = c(0, 0, 0, 0, 1, 1, 1, 1),
    time = c(0, 1, 0, 1, 0, 2, 0, 2),
    status = c(1, 0, 1, 0, 1, 0, 1, 0),
    value = c(1, 1, 1, 1, 2, 2, 2, 2)
  )
  out <- CompareAURMCs(df)
  expect_s4_class(out, "AURMC")
  expect_identical(slotNames(out), c("Arm0", "Arm1", "Contrast"))
  expect_s3_class(out@Arm0, "data.frame")
  expect_s3_class(out@Arm1, "data.frame")
  expect_s3_class(out@Contrast, "data.frame")
  expect_true("arm" %in% names(out@Arm0))
  expect_true("arm" %in% names(out@Arm1))
  expect_true(all(c("stat", "est", "lower", "upper", "p") %in% names(out@Contrast)))
  expect_equal(unique(out@Arm0$tau), 1)
  expect_equal(unique(out@Arm1$tau), 1)
})

test_that("ratio inference is unavailable for non-positive arm areas.", {
  df <- data.frame(
    idx = rep(1:4, each = 2),
    arm = rep(c(0, 0, 1, 1), each = 2),
    time = rep(c(0, 1), 4),
    status = rep(c(1, 0), 4),
    value = rep(c(-1, -1, 1, 1), each = 2)
  )
  out <- CompareAURMCs(df)
  ratio <- out@Contrast[out@Contrast$stat == "A1/A0", ]
  expect_true(is.na(ratio$se))
  expect_true(is.na(ratio$p))
})

test_that("CompareAURMCs contrast includes difference and ratio.", {
  df <- data.frame(
    idx = rep(1:4, each = 2),
    arm = rep(c(0, 0, 1, 1), each = 2),
    time = rep(c(0, 1), 4),
    status = rep(c(1, 0), 4),
    value = rep(1, 8)
  )
  out <- CompareAURMCs(df)
  stats <- unique(out@Contrast$stat)
  expect_true("A1-A0" %in% stats)
  expect_true("A1/A0" %in% stats)
})

test_that("CompareAURMCs print/show do not error.", {
  # Use differing arm values so contrast has finite se/p (avoids NaN in print)
  df <- data.frame(
    idx = c(1, 1, 2, 2, 3, 3, 4, 4),
    arm = c(0, 0, 0, 0, 1, 1, 1, 1),
    time = c(0, 1, 0, 1, 0, 2, 0, 2),
    status = c(1, 0, 1, 0, 1, 0, 1, 0),
    value = c(1, 1, 1, 1, 2, 2, 2, 2)
  )
  out <- CompareAURMCs(df)
  expect_error(print(out), NA)
  expect_error(show(out), NA)
})
