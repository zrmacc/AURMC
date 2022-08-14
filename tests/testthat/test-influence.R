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
    value = data$values
  )
  n <- est$nar[1]
  est$y <- est$nar / n
  est$d <- est$sum_value / n
  
  # Observed.
  obs <- CalcMuR(
    d = est$d,
    surv = est$surv,
    unique_times = est$time,
    y = est$y
  )
  obs <- c(obs)
  
  # Expected.
  exp <- CalcMu(
    d = est$d,
    surv = est$surv,
    unique_time = est$time,
    y = est$y
  )
  expect_equal(obs, exp)
  
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


# -----------------------------------------------------------------------------


test_that("Calculation of influence function.", {
  
  data <- data.frame(
    idx = c(1, 1, 2, 2, 3, 3, 4, 4),
    time = c(0, 1, 0, 2, 0, 3, 0, 4),
    status = c(1, 2, 1, 2, 1, 2, 1, 2),
    value = c(1, 1, 1, 1, 1, 1, 1, 1)
  )
  
  # Observed.
  influence <- InfluenceR(
    idx = data$idx,
    status = data$status,
    time = data$time,
    trunc_time = 3,
    value = data$value
  )
  
  # Expected for I1.
  est <- EstimatorR(
    idx = data$idx,
    status = data$status,
    time = data$time, 
    value = data$value,
    trunc_time = 3
  )

  n <- 4
  est$y <- est$nar / n
  est$d <- est$sum_value / n
  
  # Calculate mu.
  mu <- CalcMuR(
    d = est$d,
    surv = est$surv, 
    unique_times = est$time,
    y = est$y
  )
  mu <- as.numeric(mu)
  
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
  
  # Check I1.
  exp <- CalcI1(dm, mu, est$y)
  expect_equal(influence$i1, exp)
  
  # Expected for I2.
  value_mat <- ValueMatrixR(
    eval_times = est$time,
    idx = data$idx,
    time = data$time,
    value = data$value
  )
  exp <- CalcI2(
    d = est$d,
    surv = est$surv,
    unique_times = est$time,
    value_mat = value_mat,
    y = est$y
  )
  expect_equal(influence$i2, exp)
  
  # Expected for I3.
  risk_mat <- AtRiskMatrixR(
    eval_times = est$time,
    idx = data$idx,
    time = data$time
  )
  exp <- CalcI3(
    d = est$d,
    risk_mat = risk_mat,
    surv = est$surv,
    unique_times = est$time,
    y = est$y
  )
  expect_equal(influence$i3, exp)
  
})