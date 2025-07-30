
# Area Under the Repeated Measures Curve

<!-- badges: start -->

[![R-CMD-check](https://github.com/zrmacc/AURMC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zrmacc/AURMC/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Zachary R. McCaw <br> Updated: 2025-07-30

## Description

This package performs inference on the area under the repeated measures
curve.

## Installation

``` r
devtools::install_github(repo = 'zrmacc/AURMC')
```

## Methodological details

See the [Supplemental Proofs](docs/supplemental-proofs.pdf).

## Data

The function `GenData` simulates example data in the format expected by
this package.

``` r
library(AURMC)
set.seed(101)
data <- AURMC::GenData(
  censoring_rate = 0.5,
  death_rate = 0.5,
  n = 25,
  tau = 4
)
head(data)
```

    ##   idx      time status      value
    ## 1   1 0.0000000      1  0.5524619
    ## 2   1 1.0000000      1 -0.6749438
    ## 3   1 2.0000000      1  0.2143595
    ## 4   1 2.3638814      0  0.3107692
    ## 5   2 0.0000000      1  0.9170283
    ## 6   2 0.6309311      2 -0.2232594

The columns are:

- `idx`, the subject identifier.
- `time`, the observation time.
- `status`, coded 0 for censoring, 1 for a measurement, and 2 for a
  terminal event (e.g. death).
- `value`, the value of the measurement at the observation time.
  - The value may be continuous (as in this example) or binary.

Note that:

- All subjects must have a baseline record at `time = 0` with
  `status = 1`.
- All subjects must have an **observation-terminating** event as their
  last record, either `status = 0` for censoring or `status = 2` for a
  terminal event.
  - By default, a subject with no observation-terminating event is
    censored immediately after their last measurement.
- If a subject’s measurement value is missing (`NA`), their last value
  is carried forward.
  - If a subject’s baseline measurement is missing, it is set to zero.

## Interpolation

The function `InterpolateR` *optionally* performs linear interpolation
between measurements to enable more precise area estimation. This step
can be skipped if interpolation is not meaningful for the outcome of
interest, for example the value remains constant between observations.

``` r
# Simple example data.
data <- data.frame(
  idx = c(1, 1, 2, 2),
  status = c(1, 0, 1, 2),
  time = c(0, 1, 0, 1),
  value = c(0, 1, 0, -1)
)
show(data)
```

    ##   idx status time value
    ## 1   1      1    0     0
    ## 2   1      0    1     1
    ## 3   2      1    0     0
    ## 4   2      2    1    -1

``` r
grid <- c(0, 0.5, 1.0)
interpolated <- AURMC::Interpolate(
  data = data,
  grid = grid
)
show(interpolated)
```

    ##   idx status time value
    ## 1   1      1  0.0   0.0
    ## 2   1      1  0.5   0.5
    ## 3   1      0  1.0   1.0
    ## 4   2      1  0.0   0.0
    ## 5   2      1  0.5  -0.5
    ## 6   2      2  1.0  -1.0

## One-sample Problem

Case where the expected area is zero:

``` r
set.seed(101)
data <- AURMC::GenData(
  censoring_rate = 0.5,
  death_rate = 0.5,
  n = 25,
  tau = 4
)
AURMC::AURMC(data, tau = 2.0)
```

    ##       method tau         auc        se     lower     upper         p
    ## 1 asymptotic   2 -0.01260658 0.2577937 -0.517873 0.4926599 0.9609975

Case where the expected area is non-zero:

``` r
set.seed(101)
data <- AURMC::GenData(
  censoring_rate = 0.5,
  death_rate = 0.5,
  n = 25,
  tau = 4,
  value_mean = 1.0
)
AURMC::AURMC(data, tau = 2.0)
```

    ##       method tau    auc        se     lower    upper          p
    ## 1 asymptotic   2 1.1985 0.2735012 0.6624475 1.734553 1.1756e-05

## Two-sample Problem

Case of no true difference:

``` r
set.seed(102)

# Reference arm.
data0 <- AURMC::GenData(
  censoring_rate = 0.5,
  death_rate = 0.5,
  n = 25,
  tau = 4,
  value_mean = 1.0
)
data0$arm <- 0

# Treatment arm.
data1 <- AURMC::GenData(
  censoring_rate = 0.5,
  death_rate = 0.5,
  n = 25,
  tau = 4,
  value_mean = 1.0
)
data1$arm <- 1
data1$idx <- nrow(data0) + data1$idx

# Overall data set.
data <- rbind(data0, data1)

AURMC::CompareAURMCs(data)
```

    ## Arm 0:
    ##       method arm tau auc    se lower upper        p
    ## 1 asymptotic   0 3.6   1 0.227 0.556  1.45 1.01e-05
    ## 
    ## 
    ## Arm 1:
    ##       method arm tau  auc    se lower upper        p
    ## 1 asymptotic   1 3.6 1.34 0.278 0.793  1.88 1.47e-06
    ## 
    ## 
    ## Contrast:
    ##       method  stat   est    se  lower upper     p
    ## 1 asymptotic A1-A0 0.336 0.358 -0.366  1.04 0.348
    ## 2 asymptotic A1/A0 1.340 0.411  0.732  2.44 0.346

Case of a true differnece:

``` r
set.seed(102)

# Reference arm.
data0 <- AURMC::GenData(
  censoring_rate = 0.5,
  death_rate = 0.5,
  n = 25,
  tau = 4,
  value_mean = 1.0
)
data0$arm <- 0

# Treatment arm.
data1 <- AURMC::GenData(
  censoring_rate = 0.5,
  death_rate = 0.5,
  n = 25,
  tau = 4,
  value_mean = 2.0
)
data1$arm <- 1
data1$idx <- nrow(data0) + data1$idx

# Overall data set.
data <- rbind(data0, data1)

AURMC::CompareAURMCs(data)
```

    ## Arm 0:
    ##       method arm tau auc    se lower upper        p
    ## 1 asymptotic   0 3.6   1 0.227 0.556  1.45 1.01e-05
    ## 
    ## 
    ## Arm 1:
    ##       method arm tau  auc    se lower upper        p
    ## 1 asymptotic   1 3.6 2.53 0.458  1.63  3.43 3.44e-08
    ## 
    ## 
    ## Contrast:
    ##       method  stat  est    se lower upper       p
    ## 1 asymptotic A1-A0 1.53 0.511 0.526  2.53 0.00281
    ## 2 asymptotic A1/A0 2.53 0.733 1.430  4.46 0.00140
