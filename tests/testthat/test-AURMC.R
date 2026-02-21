# Main AURMC() tests (smoke, structure, bootstrap)

test_that("AURMC returns expected structure.", {
  df <- data.frame(
    idx = c(1, 1, 2, 2, 2),
    time = c(0, 1, 0, 1, 2),
    status = c(1, 0, 1, 1, 0),
    value = c(1, 1, 1, 1, 1)
  )
  out <- AURMC(df)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("method", "tau", "auc", "se", "lower", "upper", "p"))
  expect_equal(out$method, "asymptotic")
  expect_equal(nrow(out), 1L)
  expect_true(is.numeric(out$auc) && length(out$auc) == 1L)
  expect_true(out$lower <= out$auc && out$auc <= out$upper)
})

test_that("AURMC with perturbations returns two rows.", {
  df <- data.frame(
    idx = c(1, 1, 2, 2, 2),
    time = c(0, 1, 0, 1, 2),
    status = c(1, 0, 1, 1, 0),
    value = c(1, 1, 1, 1, 1)
  )
  out <- AURMC(df, perturbations = 50L, random_state = 1L)
  expect_equal(nrow(out), 2L)
  expect_equal(out$method, c("asymptotic", "bootstrap"))
  expect_true(all(is.finite(out$auc)))
  expect_true(all(is.finite(out$se)))
})

test_that("AURMC respects tau.", {
  df <- data.frame(
    idx = c(1, 1, 1, 1),
    time = c(0, 1, 2, 3),
    status = c(1, 1, 1, 0),
    value = c(1, 1, 1, 1)
  )
  out_full <- AURMC(df, tau = 3)
  out_trunc <- AURMC(df, tau = 1.5)
  expect_equal(out_full$tau, 3)
  expect_equal(out_trunc$tau, 1.5)
  expect_true(out_trunc$auc < out_full$auc)
})

test_that("AURMC with custom column names.", {
  df <- data.frame(
    id = c(1, 1, 2, 2),
    t = c(0, 1, 0, 2),
    s = c(1, 0, 1, 0),
    v = c(1, 1, 1, 1)
  )
  out <- AURMC(df, idx_name = "id", time_name = "t", status_name = "s", value_name = "v")
  expect_equal(nrow(out), 1L)
  expect_true(is.finite(out$auc))
})
