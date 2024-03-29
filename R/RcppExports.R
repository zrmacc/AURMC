# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Tabulate Value Matrix R
#'
#' Tabulate \eqn{D_{i}(t)} with subjects as rows and evaluation times as columns.
#'  
#' @param eval_times Evaluation times.
#' @param idx Unique subject index. 
#' @param time Observation time.
#' @param value Observation value.
#' @return Numeric matrix.
ValueMatrixR <- function(eval_times, idx, time, value) {
    .Call(`_AURMC_ValueMatrixR`, eval_times, idx, time, value)
}

#' Tabulate At-Risk Matrix R
#'
#' Construct a subject (row) by evaluation time (col) matrix.
#'  
#' @param eval_times Evaluation times.
#' @param idx Unique subject index. 
#' @param time Observation time.
#' @return Numeric matrix.
AtRiskMatrixR <- function(eval_times, idx, time) {
    .Call(`_AURMC_AtRiskMatrixR`, eval_times, idx, time)
}

#' Tabulate Kaplan Meier R
#'
#' Constructs a matrix with evaluation times as rows, and 4 columns:
#' \itemize{
#' \item{time}{Evaluation times.}
#' \item{nar}{Number at risk.}
#' \item{surv}{Survival probability.}
#' \item{haz}{Hazard.}
#' }
#'  
#' @param eval_times Evaluation times.
#' @param idx Unique subject index.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
#' @param time Observation time.
#' @return Data.frame.
KaplanMeierR <- function(eval_times, idx, status, time) {
    .Call(`_AURMC_KaplanMeierR`, eval_times, idx, status, time)
}

#' Tabulate Estimator R
#'  
#' @param idx Unique subject index. 
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
#' @param time Observation time.
#' @param value Observation value.
#' @param eval_times Evalulation times. If omitted, defaults to the
#' unique values of time.
#' @param replace_na Replace NaN with zero? Default: FALSE.
#' @param return_auc Return the AUC? Default: FALSE.
#' @param trunc_time Truncation time? Optional. If omitted, defaults
#' to the maximum evaluation time.
#' @return Data.frame.
EstimatorR <- function(idx, status, time, value, eval_times = NULL, replace_na = FALSE, return_auc = FALSE, trunc_time = NULL) {
    .Call(`_AURMC_EstimatorR`, idx, status, time, value, eval_times, replace_na, return_auc, trunc_time)
}

#' Draw Bootstrap R
#'  
#' @param idx Unique subject index. 
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
#' @param time Observation time.
#' @param value Observation value.
#' @return Numeric matrix.
DrawBootstrapR <- function(idx, status, time, value) {
    .Call(`_AURMC_DrawBootstrapR`, idx, status, time, value)
}

#' Bootstrap Samples R
#'
#' Returns a bootstrap (row) by statistic (col) matrix. If \code{return_auc = TRUE},
#' the output has a single column, the AUC. If \code{return_auc = FALSE}, the output
#' has a single column for each evaluation time.
#'  
#' @param boot Bootstrap replicates.
#' @param eval_times Evaluation times.
#' @param idx Unique subject index. 
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
#' @param time Observation time.
#' @param value Observation value.
#' @param replace_na Replace NaN with zero?
#' @param return_auc Return the AUC?
#' @param trunc_time Truncation time? Optional. If omitted, defaults
#' to the maximum evaluation time.
#' @return Numeric matrix.
BootstrapSamplesR <- function(boot, eval_times, idx, status, time, value, replace_na = FALSE, return_auc = FALSE, trunc_time = NULL) {
    .Call(`_AURMC_BootstrapSamplesR`, boot, eval_times, idx, status, time, value, replace_na, return_auc, trunc_time)
}

#' Calculate Mu R
#'
#' Evaluate \eqn{\mu(t; \tau) = \int_{t}^{\tau}{S(u)d(u)/y(u)}du}.
#' 
#' @param d Value of d(t) at each time point.
#' @param surv Value of S(t) at each time point.
#' @param unique_times Unique values of time t.
#' @param y Value of y(t) at each time point.
#' @return Numeric vector of \eqn{\mu(t; tau)}.
CalcMuR <- function(d, surv, unique_times, y) {
    .Call(`_AURMC_CalcMuR`, d, surv, unique_times, y)
}

#' Calculate Martingale
#'
#' Calculate \eqn{dM_{i}(t) = dN_{i}(t) - Y_{i}(t)d\Lambda(t)}.
#' 
#' @param haz Value of the hazard at each unique time.
#' @param idx Subject index.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
#' @param time Subject observation times.
#' @param unique_times Unique times at which to obtain the martingale.
#' @return Matrix with subjects as rows and unique times as columns.
CalcMartingaleR <- function(haz, idx, status, time, unique_times) {
    .Call(`_AURMC_CalcMartingaleR`, haz, idx, status, time, unique_times)
}

#' Influence Function R
#' 
#' Influence function contributions for the AUC. Includes the three component
#' integrals (`i1`, `i2`, `i3`) and the overall influence `psi`.
#'  
#' @param idx Unique subject index. 
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
#' @param time Observation time.
#' @param trunc_time Truncation time? Optional. If omitted, defaults
#' to the maximum evaluation time.
#' @param value Observation value.
#' @return Data.frame.
InfluenceR <- function(idx, status, time, trunc_time, value) {
    .Call(`_AURMC_InfluenceR`, idx, status, time, trunc_time, value)
}

#' Perturbation R
#'  
#' Generates realizations of \eqn{\frac{1}{n}\sum_{i=1}^{n}\psi_{i}w_{i}},
#' where \eqn{\psi_{i}} is the influence function for the ith subject and the
#' \eqn{w_{i}} are IID random weights.
#' 
#' @section Notes:
#' The random seed should be set in R prior to calling this function.
#'  
#' @param idx Unique subject index. 
#' @param perturbations Number of perturbations
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
#' @param time Observation time.
#' @param trunc_time Truncation time? Optional. If omitted, defaults
#' to the maximum evaluation time.
#' @param value Observation value.
#' @return Numeric vector.
PerturbationR <- function(idx, perturbations, status, time, trunc_time, value) {
    .Call(`_AURMC_PerturbationR`, idx, perturbations, status, time, trunc_time, value)
}

#' Interpolation R
#'  
#' Linearly interpolations between each subject's measurements.
#' The input data should contain no missing values. 
#'  
#' @param grid Grid of unique points at which to interpolate.
#' @param idx Unique subject index. 
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
#' @param time Observation time.
#' @param value Observation value.
#' @return Data.frame.
InterpolateR <- function(grid, idx, status, time, value) {
    .Call(`_AURMC_InterpolateR`, grid, idx, status, time, value)
}

