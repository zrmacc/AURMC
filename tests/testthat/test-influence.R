test_that("Calculation of mu.", {
  
  withr::local_seed(101)
  data <- AURMC::GenData(
    censoring_rate = 0.50,
    death_rate = 0.25,
    n = 10, 
    tau = 5
  )
  
  est <- EstimatorR(
    idx = data$idx,
    status = data$status,
    time = data$time,
    trunc_time = 1.0,
    value = data$value
  )
  
  for (int_method in c("left", "right", "trapezoid")) {
    obs <- CalcMuR(
      d = est$d,
      surv = est$surv,
      unique_times = est$time,
      y = est$y,
      int_method = int_method
    )
    exp <- CalcMu(
      d = est$d,
      surv = est$surv,
      unique_times = est$time,
      y = est$y,
      int_method = int_method
    )
    expect_equal(c(obs), exp)
  }
  
})

test_that("Influence tail excludes the contemporaneous hazard jump.", {
  unique_times <- c(0, 1, 2)
  d <- surv <- y <- rep(1, 3)

  expected <- list(
    left = c(1, 0, 0),
    right = c(2, 1, 0),
    trapezoid = c(1.5, 0.5, 0)
  )

  for (int_method in names(expected)) {
    expect_equal(
      as.numeric(CalcMuR(d, surv, unique_times, y, int_method)),
      expected[[int_method]]
    )
    expect_equal(
      CalcMu(d, surv, unique_times, y, int_method),
      expected[[int_method]]
    )
  }
})


# -----------------------------------------------------------------------------


test_that("Calculation of martingales.", {
  
  # Function to calculate Kaplan-Meier curve.
  GetKM <- function(data) {
    KaplanMeierR(
      eval_times = sort(unique(data$time)),
      idx = data$idx,
      status = data$status,
      time = data$time
    )
  }
  
  # Case 1: no deaths.
  data <- data.frame(
    idx = c(1, 1, 2, 2),
    time = c(0, 1, 0, 2),
    status = c(1, 0, 1, 0)
  )
  km <- GetKM(data)
  
  # Observed.
  obs <- CalcMartingaleR(
    haz = km$haz,
    idx = data$idx,
    time = data$time,
    status = data$status,
    unique_times = c(0, 1, 2)
  )
  
  # Expected.
  exp <- array(0, dim = c(2, 3))
  expect_equal(obs, exp)
  
  # Case 2.
  data <- data.frame(
    idx = c(1, 1, 2, 2),
    time = c(0, 1, 0, 2),
    status = c(1, 2, 1, 0)
  )
  km <- GetKM(data)
  
  # Observed.
  obs <- CalcMartingaleR(
    haz = km$haz,
    idx = data$idx,
    time = data$time,
    status = data$status,
    unique_times = c(0, 1, 2)
  )
  
  # Expected.
  exp <- array(0, dim = c(2, 3))
  exp[1, 2] <- 1.0 - 0.5
  exp[2, 2] <- 0.0 - 0.5
  expect_equal(obs, exp)
  
  # Case 3.
  data <- data.frame(
    idx = c(1, 1, 2, 2, 3, 3),
    time = c(0, 1, 0, 2, 0, 3),
    status = c(1, 0, 1, 2, 1, 0)
  )
  km <- GetKM(data)
  
  # Observed.
  obs <- CalcMartingaleR(
    haz = km$haz,
    idx = data$idx,
    time = data$time,
    status = data$status,
    unique_times = c(0, 1, 2, 3)
  )
  
  # Expected.
  exp <- array(0, dim = c(3, 4))
  exp[1, 3] <- 0.0 - 0.0
  exp[2, 3] <- 1.0 - 0.5
  exp[3, 3] <- 0.0 - 0.5
  expect_equal(obs, exp)  
  
})

test_that("A censoring record tied with a measurement is not treated as death.", {
  data <- data.frame(
    idx = c(1, 1, 1, 2, 2),
    time = c(0, 1, 1, 0, 2),
    status = c(1, 0, 1, 1, 0)
  )
  km <- KaplanMeierR(
    eval_times = c(0, 1, 2),
    idx = data$idx,
    status = data$status,
    time = data$time
  )
  observed <- CalcMartingaleR(
    haz = km$haz,
    idx = data$idx,
    status = data$status,
    time = data$time,
    unique_times = c(0, 1, 2)
  )
  expect_equal(observed, matrix(0, nrow = 2, ncol = 3))
})


# -----------------------------------------------------------------------------


test_that("Calculation of influence function respects integration method.", {
  
  data <- data.frame(
    idx = c(1, 1, 2, 2, 3, 3, 4, 4),
    time = c(0, 1, 0, 2, 0, 3, 0, 4),
    status = c(1, 2, 1, 2, 1, 2, 1, 2),
    value = c(1, 1, 1, 1, 1, 1, 1, 1)
  )
  
  # Expected for I1.
  est <- EstimatorR(
    idx = data$idx,
    status = data$status,
    time = data$time, 
    value = data$value,
    trunc_time = 3
  )
  
  # Calculate Kaplan-Meier.
  km <- KaplanMeierR(
    eval_times = est$time,
    idx = data$idx,
    status = data$status,
    time = data$time
  )
  
  # Calculate martingales.
  dm <- CalcMartingaleR(
    haz = km$haz,
    idx = data$idx,
    status = data$status,
    time = data$time,
    unique_times = est$time
  )
  
  value_mat <- ValueMatrixR(
    eval_times = est$time,
    idx = data$idx,
    time = data$time,
    value = data$value
  )
  risk_mat <- AtRiskMatrixR(
    eval_times = est$time,
    idx = data$idx,
    time = data$time
  )
  for (int_method in c("left", "right", "trapezoid")) {
    influence <- InfluenceR(
      idx = data$idx,
      status = data$status,
      time = data$time,
      trunc_time = 3,
      value = data$value,
      int_method = int_method
    )

    mu <- CalcMuR(
      d = est$d,
      surv = est$surv,
      unique_times = est$time,
      y = est$y,
      int_method = int_method
    )
    expect_equal(influence$i1, CalcI1(dm, as.numeric(mu), est$y))
    expect_equal(
      influence$i2,
      CalcI2(est$d, est$surv, est$time, value_mat, est$y, int_method)
    )
    expect_equal(
      influence$i3,
      CalcI3(est$d, risk_mat, est$surv, est$time, est$y, int_method)
    )
  }
  
})
