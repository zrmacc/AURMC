# Area Under the Repeated Measures Curve

Zachary R. McCaw <br>
Updated: 2022-08-14



## Description

This package performs inference on the area under the repeated measures curve.

## Installation


```r
devtools::install_github(repo = 'zrmacc/AURMC')
```

## Data

The function `GenData` simulates example data in the format expected by this package. 


```r
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

```
##   idx      time status      value
## 1   1 0.0000000      1  0.5524619
## 2   1 1.0000000      1 -0.6749438
## 3   1 2.0000000      1  0.2143595
## 4   1 2.3638814      0  0.3107692
## 5   2 0.0000000      1  0.9170283
## 6   2 0.6309311      2 -0.2232594
```

The columns are:

* `idx`, the subject identifier. 
* `time`, the observation time. 
* `status`, coded 0 for censoring, 1 for a measurement, and 2 for a terminal event (e.g. death).
* `value`, the value of the measurement at the observation time.

Note that:

* All subjects must have a baseline record at `time = 0` with `status = 1`.
* All subjects must have an *observation-terminating* event as their last record, either `status = 0` for censoring or `status = 2` for a terminal event.
  - By default, a subject with no observation-terminating event is censored immediately after their last measurement.
* If a subject's measurement value is missing (`NA`), their last value is carried forward.
  - If a subject's baseline measurement is missing, it is set to zero.

## Usage

Case where the expected area is zero:

```r
set.seed(101)
data <- AURMC::GenData(
  censoring_rate = 0.5,
  death_rate = 0.5,
  n = 25,
  tau = 4
)
AURMC::AURMC(data, tau = 2.0)
```

```
##   tau         auc        se     lower     upper         p
## 1   2 -0.01260658 0.2577937 -0.517873 0.4926599 0.9609975
```

Case where the expected area is non-zero:

```r
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

```
##   tau    auc        se     lower    upper          p
## 1   2 1.1985 0.2735012 0.6624475 1.734553 1.1756e-05
```
