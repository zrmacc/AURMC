% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TwoSample.R
\name{CompareAURMCs}
\alias{CompareAURMCs}
\title{Compare Areas Under Repeated Measures Curves}
\usage{
CompareAURMCs(
  data,
  alpha = 0.05,
  arm_name = "arm",
  censor_after_last = TRUE,
  idx_name = "idx",
  perturbations = NULL,
  random_state = 0,
  status_name = "status",
  tau = NULL,
  time_name = "time",
  value_name = "value"
)
}
\arguments{
\item{data}{Data.frame.}

\item{alpha}{Type I error.}

\item{arm_name}{Name of the column containing treatment arm. Coded as 0 for
reference, 1 for treatment.}

\item{censor_after_last}{Introduce censoring after the last event \emph{if} no
observation-terminating event is present.}

\item{idx_name}{Name of column containing a unique subject index.}

\item{perturbations}{Number of perturbations to use for bootstrap inference.
If \code{NULL}, only analytical inference is performed.}

\item{random_state}{Seed to ensure perturbations are reproducible.}

\item{status_name}{Name of column containing the status. Must be coded as 0
for censoring, 1 for a measurement, 2 for death. Each subject should have
an observation-terminating event, either censoring or death.}

\item{tau}{Truncation time.}

\item{time_name}{Name of column containing the observation time.}

\item{value_name}{Name of the column containing the measurement.}
}
\value{
Data.frame.
}
\description{
Compare Areas Under Repeated Measures Curves
}
