// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------
// Utilities.
// ----------------------------------------------------------------------------

// Check if Value is in Vector
// 
// @param a Value to search for.
// @param b Vector to search.
// @return bool.

bool IsIn(const double &a, const arma::vec &b) {
  
  for(int i=0; i<b.size(); i++) {
    if(b(i) == a) {
      return true;
    }
  }
  return false;
}


// Union
// 
// @param a Vector.
// @param b Vector.
// @return Vector.
arma::colvec Union(const arma::colvec &a, arma::colvec b) {
  
  // Sets all elements of b that are in a to a value that is known
  // to be in a. Then concatenates a and b and filters to unique values.
  
  double a0 = a(0);
  for(int i=0; i<b.size(); i++) {
    if(IsIn(b(i), a)){
      b(i) = a0;
    }
  }

  arma::colvec out = arma::unique(arma::join_cols(a, b));
  return out;
}


// Truncate
//
// @param time Vector of time points.
// @param tau Truncation time.
// @return Truncated vector.
arma::colvec Truncate(const arma::colvec &time, const double tau) {
  arma::colvec unique_times = arma::unique(time);
  for(int i=0; i<unique_times.size(); i++){
    if(unique_times(i) > tau){
      unique_times(i) = tau;
    }
  }
  arma::colvec out = arma::unique(unique_times);
  return out;
}


// ----------------------------------------------------------------------------
// Value matrix.
// ----------------------------------------------------------------------------

//' Tabulate Value Matrix R
//'
//' Tabulate \eqn{D_{i}(t)} with subjects as rows and evaluation times as columns.
//'  
//' @param eval_times Evaluation times.
//' @param idx Unique subject index. 
//' @param time Observation time.
//' @param value Observation value.
//' @return Numeric matrix.
// [[Rcpp::export]]

SEXP ValueMatrixR(
    const arma::colvec eval_times,
    const arma::colvec idx,
    const arma::colvec time,
    const arma::colvec value
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Create a subject by evaluation times matrix, where
  // D(i, t) is subject i's current value at time t.
  const int n_times = eval_times.size();
  arma::mat dmat = arma::zeros(n, n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Current subject.
    arma::colvec subj_times = time.elem(arma::find(idx == unique_idx(i)));
    arma::colvec subj_values = value.elem(arma::find(idx == unique_idx(i)));
    
    // Initialize.
    double current_value = 0;
    double next_value = 0;

    // Loop over evaluation times.
    for(int j=0; j<n_times; j++) {

      double current_time = eval_times(j);

      if (arma::any(subj_times <= current_time)) {
        next_value = arma::as_scalar(
          subj_values.elem(arma::find(subj_times <= current_time, 1, "last"))
        );
      }

      if (not std::isnan(next_value)) {
        current_value = next_value;
      }

      dmat(i,j) = current_value;

      if (arma::all(subj_times <= current_time)) {
        break;
      }
    }
  }
  return Rcpp::wrap(dmat);
}


// ----------------------------------------------------------------------------


// Tabulate Value Matrix Cpp
//
// Tabulate \eqn{D_{i}(t)} with subjects as rows and evaluation times as columns.
//  
// @param eval_times Evaluation times.
// @param idx Unique subject index. 
// @param time Observation time.
// @param value Observation value.
// @return Numeric matrix.

arma::mat ValueMatrixCpp(
    const arma::colvec eval_times,
    const arma::colvec idx,
    const arma::colvec time,
    const arma::colvec value
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Create a subject by evaluation times matrix, where
  // D(i, t) is subject i's current value at time t.
  const int n_times = eval_times.size();
  arma::mat dmat = arma::zeros(n, n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Current subject.
    arma::colvec subj_times = time.elem(arma::find(idx == unique_idx(i)));
    arma::colvec subj_values = value.elem(arma::find(idx == unique_idx(i)));
    
    // Initialize.
    double current_value = 0;
    double next_value = 0;

    // Loop over evaluation times.
    for(int j=0; j<n_times; j++) {

      double current_time = eval_times(j);

      if (arma::any(subj_times <= current_time)) {
        next_value = arma::as_scalar(
          subj_values.elem(arma::find(subj_times <= current_time, 1, "last"))
        );
      }

      if (not std::isnan(next_value)) {
        current_value = next_value;
      }

      dmat(i,j) = current_value;

      if (arma::all(subj_times <= current_time)) {
        break;
      }
    }
  }
  return dmat;
}


// ----------------------------------------------------------------------------
// At risk matrix.
// ----------------------------------------------------------------------------


//' Tabulate At-Risk Matrix R
//'
//' Construct a subject (row) by evaluation time (col) matrix.
//'  
//' @param eval_times Evaluation times.
//' @param idx Unique subject index. 
//' @param time Observation time.
//' @return Numeric matrix.
// [[Rcpp::export]]

SEXP AtRiskMatrixR(
    const arma::colvec eval_times,
    const arma::colvec idx,
    const arma::colvec time
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Create a subject by evaluation times matrix, where
  // Y(i, t) is 1 if subject i is at-risk at time t.
  const int n_times = eval_times.size();
  arma::mat y = arma::zeros(n, n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Current subject.
    arma::colvec subj_times = time.elem(arma::find(idx == unique_idx(i)));
    
    // Subject is at risk until time > the subject's last time.
    double subj_last_time = subj_times.max();

    // Loop over times.
    for(int j=0; j<n_times; j++) {
      
      double current_time = eval_times(j);

      if(current_time <= subj_last_time) {
        y(i, j) = 1;
      } else {
        break;
      }
    }
  }
  return Rcpp::wrap(y);
}


// ----------------------------------------------------------------------------


// Tabulate At-Risk Matrix Cpp
//
// Construct a subject (row) by evaluation time (col) matrix.
//  
// @param eval_times Evaluation times.
// @param idx Unique subject index. 
// @param time Observation time.
// @return Numeric matrix.

arma::mat AtRiskMatrixCpp(
    const arma::colvec eval_times,
    const arma::colvec idx,
    const arma::colvec time
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Create a subject by evaluation times matrix, where
  // Y[i, t] is 1 if subject i is at-risk at time t.
  const int n_times = eval_times.size();
  arma::mat y = arma::zeros(n, n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Current subject.
    arma::colvec subj_times = time.elem(arma::find(idx == unique_idx(i)));
    
    // Subject is at risk until time > the subject's last time.
    double subj_last_time = subj_times.max();
    
    // Loop over times.
    for(int j=0; j<n_times; j++) {
      
      double current_time = eval_times(j);

      if(current_time <= subj_last_time) {
        y(i, j) = 1;
      } else {
        break;
      }
    }
  }
  return y;
}


// ----------------------------------------------------------------------------
// Kaplan Meier.
// ----------------------------------------------------------------------------


//' Tabulate Kaplan Meier R
//'
//' Constructs a matrix with evaluation times as rows, and 4 columns:
//' \itemize{
//' \item{time}{Evaluation times.}
//' \item{nar}{Number at risk.}
//' \item{surv}{Survival probability.}
//' \item{haz}{Hazard.}
//' }
//'  
//' @param eval_times Evaluation times.
//' @param idx Unique subject index.
//' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
//' @param time Observation time.
//' @return Data.frame.
// [[Rcpp::export]]

SEXP KaplanMeierR(
    const arma::colvec eval_times,
    const arma::colvec idx,
    const arma::colvec status,
    const arma::colvec time
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Unique times.
  arma::colvec unique_times = Union(eval_times, arma::unique(time));
  const int n_unique_time = unique_times.size();
  
  // Censoring, death, and at risk counts.
  arma::colvec censor(n_unique_time);
  arma::colvec death(n_unique_time);
  arma::colvec nar(n_unique_time);
  
  // Loop over unique times.
  double current_nar = n;

  for(int i=0; i<n_unique_time; i++) {
    
    double current_time = unique_times(i);
    const arma::colvec current_status = status.elem(arma::find(time == current_time));
    
    nar(i) = current_nar;
    censor(i) = arma::sum(current_status == 0.0);
    death(i) = arma::sum(current_status == 2.0);
    
    // Update NAR.
    current_nar -= censor(i) + death(i);
  }
  
  // Hazard.
  const arma::colvec haz = death / nar;
  
  // Survival probability.
  arma::colvec surv = arma::cumprod(1 - haz);
  
  // Restrict to evaluation times.
  const int n_eval_time = eval_times.size();
  arma::colvec nar_out(n_eval_time);
  arma::colvec haz_out(n_eval_time);
  arma::colvec surv_out(n_eval_time);

  int pointer = 0;
  for(int i=0; i<n_unique_time; i++) {
    
    double current_time = unique_times(i);
    if(IsIn(current_time, eval_times)) {
      nar_out(pointer) = nar(i);
      haz_out(pointer) = haz(i);
      surv_out(pointer) = surv(i);
      pointer += 1;
    }
    
  }

  // Output.
  return Rcpp::DataFrame::create(
    Rcpp::Named("time")=eval_times,
    Rcpp::Named("nar")=nar_out,
    Rcpp::Named("surv")=surv_out,
    Rcpp::Named("haz")=haz_out
  );
}


// ----------------------------------------------------------------------------


// Tabulate Kaplan Meier Cpp
//  
// Constructs a matrix with evaluation times as rows, and 4 columns:
// \itemize{
// \item{time}{Evaluation times.}
// \item{nar}{Number at risk.}
// \item{surv}{Survival probability.}
// \item{haz}{Hazard.}
// }
//
// @param eval_times Evaluation times.
// @param idx Unique subject index.
// @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
// @param time Observation time.
// @return Numeric matrix.

arma::mat KaplanMeierCpp(
    const arma::colvec eval_times,
    const arma::colvec idx,
    const arma::colvec status,
    const arma::colvec time
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Unique times.
  arma::colvec unique_times = Union(eval_times, arma::unique(time));
  const int n_unique_time = unique_times.size();
  
  // Censoring, death, and at risk counts.
  arma::colvec censor(n_unique_time);
  arma::colvec death(n_unique_time);
  arma::colvec nar(n_unique_time);
  
  // Loop over unique times.
  double current_nar = n;

  for(int i=0; i<n_unique_time; i++) {
    
    double current_time = unique_times(i);
    const arma::colvec current_status = status.elem(arma::find(time == current_time));
    
    nar(i) = current_nar;
    censor(i) = arma::sum(current_status == 0.0);
    death(i) = arma::sum(current_status == 2.0);
    
    // Update NAR.
    current_nar -= censor(i) + death(i);
  }
  
  // Hazard.
  const arma::colvec haz = death / nar;
  
  // Survival probability.
  arma::colvec surv = arma::cumprod(1 - haz);
  
  // Restrict to evaluation times.
  const int n_eval_time = eval_times.size();
  arma::colvec nar_out(n_eval_time);
  arma::colvec haz_out(n_eval_time);
  arma::colvec surv_out(n_eval_time);

  int pointer = 0;
  for(int i=0; i<n_unique_time; i++) {
    
    double current_time = unique_times(i);
    if(IsIn(current_time, eval_times)) {
      nar_out(pointer) = nar(i);
      haz_out(pointer) = haz(i);
      surv_out(pointer) = surv(i);
      pointer += 1;
    }
    
  }
  
  // Output.
  arma::mat out = arma::join_rows(eval_times, nar_out, surv_out, haz_out);
  return out;
}


// ----------------------------------------------------------------------------
// Estimators.
// ----------------------------------------------------------------------------

//' Tabulate Estimator R
//'  
//' @param idx Unique subject index. 
//' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
//' @param time Observation time.
//' @param value Observation value.
//' @param eval_times Evalulation times. If omitted, defaults to the
//' unique values of time.
//' @param replace_na Replace NaN with zero? Default: FALSE.
//' @param return_auc Return the AUC? Default: FALSE.
//' @param trunc_time Truncation time? Optional. If omitted, defaults
//' to the maximum evaluation time.
//' @return Data.frame.
// [[Rcpp::export]]

SEXP EstimatorR(
  const arma::colvec idx,
  const arma::colvec status,
  const arma::colvec time,
  const arma::colvec value,
  const Rcpp::Nullable<Rcpp::NumericVector> eval_times=R_NilValue,
  const bool replace_na=false,
  const bool return_auc=false,
  const Rcpp::Nullable<double> trunc_time=R_NilValue
){

  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();

  // Unique times.
  // Armadillo's unique function sorts the values in ascending order.
  arma::colvec unique_times;
  if(eval_times.isNotNull()) {
    unique_times = Rcpp::as<arma::colvec>(eval_times);
    unique_times = arma::unique(unique_times);
  } else {
    unique_times = arma::unique(time);
  }

  // Truncate.
  if(trunc_time.isNotNull()) {
    double tau = Rcpp::as<double>(trunc_time);
    unique_times = Truncate(unique_times, tau);
  } 
  const int n_times = unique_times.size();

  // Tabulate D_{i}(t).
  arma::mat value_mat = ValueMatrixCpp(unique_times, idx, time, value);

  // Estimate d(t).
  arma::colvec d = arma::trans(arma::sum(value_mat, 0)) / n;
  
  // Kaplan-Meier.
  arma::mat km_mat = KaplanMeierCpp(unique_times, idx, status, time);
  if(replace_na) {
  	km_mat.replace(arma::datum::nan, 0);
  }

  const arma::colvec nar = km_mat.col(1);
  const arma::colvec surv = km_mat.col(2);
  const arma::colvec haz = km_mat.col(3);

  // Proportion at risk.
  const arma::colvec y = nar / n;

  // Expectation E{D(t)}.
  arma::colvec exp = surv % (d / y);
  if(replace_na) {
  	exp.replace(arma::datum::nan, 0);
  }

  if(return_auc) {

    const arma::colvec delta_t = arma::diff(unique_times);
    const arma::colvec integrand = exp.subvec(0, n_times - 2);
    double auc = arma::sum(integrand % delta_t);
    return Rcpp::wrap(auc);

  }

  // Output.
  return Rcpp::DataFrame::create(
    Rcpp::Named("time")=unique_times,
    Rcpp::Named("nar")=nar,
    Rcpp::Named("y")=y,
    Rcpp::Named("haz")=haz,
    Rcpp::Named("surv")=surv,
    Rcpp::Named("d")=d,
    Rcpp::Named("exp")=exp
  );
}


// ----------------------------------------------------------------------------


// Tabulate Estimator Cpp
//
// If AUC is requested, returns 1x1 AUC. Otherwise, returns the estimator
// at each evaluation time.
//  
// @param eval_times Evaluation times.
// @param idx Unique subject index. 
// @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
// @param time Observation time.
// @param value Observation value.
// @param replace_na Replace NaN with zero?
// @param return_auc Return the AUC?
// @return Numeric vector.

arma::mat EstimatorCpp(
  arma::colvec eval_times,
  const arma::colvec idx,
  const arma::colvec status,
  const arma::colvec time,
  const double trunc_time,
  const arma::colvec value,
  const bool replace_na=false,
  const bool return_auc=false
){

  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Evaluation times.
  eval_times = Truncate(eval_times, trunc_time);
  const int n_times = eval_times.size();

  // Tabulate D_{i}(t).
  arma::mat value_mat = ValueMatrixCpp(eval_times, idx, time, value);

  // Estimate d(t).
  arma::colvec d = arma::trans(arma::sum(value_mat, 0)) / n;
  
  // Construct Kaplan-Meier curve: evalulation time (row) by 2.
  arma::mat km_mat =  KaplanMeierCpp(eval_times, idx, status, time);
  if(replace_na) {
  	km_mat.replace(arma::datum::nan, 0);
  }
  arma::colvec y = km_mat.col(1) / n;
  arma::colvec surv = km_mat.col(2);
  arma::colvec haz = km_mat.col(3);

  // Expectation: E{D(t)}.
  arma::colvec exp = surv % (d / y);
  if(replace_na) {
  	exp.replace(arma::datum::nan, 0);
  }

  if(return_auc) {

    const arma::colvec delta_t = arma::diff(eval_times);
    const arma::colvec integrand = exp.subvec(0, n_times - 2);
    double auc = arma::sum(integrand % delta_t);
    arma::mat out(1, 1);
    out(0, 0) = auc;
    return out;

  } else {

    arma::mat out = arma::join_rows(eval_times, y, haz, surv);
    out = arma::join_rows(out, d, exp);
  	return out;

  }
}


// ----------------------------------------------------------------------------
// Bootstrap
// ----------------------------------------------------------------------------

//' Draw Bootstrap R
//'  
//' @param idx Unique subject index. 
//' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
//' @param time Observation time.
//' @param value Observation value.
//' @return Numeric matrix.
// [[Rcpp::export]]

SEXP DrawBootstrapR(
  const arma::colvec idx,
  const arma::colvec status,
  const arma::colvec time,
  const arma::colvec value
){

  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();

  // Bind columns to form matrix.
  arma::mat data = arma::join_rows(idx, status, time, value);

  // Bootstrap.
  arma::mat boot(0, 4);
  for(int i=0; i<n; i++) {

    const int draw = arma::randi(arma::distr_param(0, n - 1));
    const int draw_idx = unique_idx(draw);
    arma::mat slice = data.rows(arma::find(idx == draw_idx));
    slice.col(0) = i * arma::ones(slice.n_rows);
    boot = arma::join_cols(boot, slice);

  }

  return Rcpp::wrap(boot);
}


// ----------------------------------------------------------------------------


// Draw Bootstrap Cpp
//  
// @param idx Unique subject index. 
// @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
// @param time Observation time.
// @param value Observation value.
// @return Numeric matrix.

arma::mat DrawBootstrapCpp(
  const arma::colvec idx,
  const arma::colvec status,
  const arma::colvec time,
  const arma::colvec value
){

  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();

  // Bind columns to form matrix.
  arma::mat data = arma::join_rows(idx, status, time, value);
  

  // Bootstrap.
  arma::mat boot(0, 4);
  for(int i=0; i<n; i++) {

    const int draw = arma::randi(arma::distr_param(0, n - 1));
    const int draw_idx = unique_idx(draw);
    arma::mat slice = data.rows(arma::find(idx == draw_idx));
    slice.col(0) = i * arma::ones(slice.n_rows);
    boot = arma::join_cols(boot, slice);

  }

  return boot;
}


// ----------------------------------------------------------------------------


//' Bootstrap Samples R
//'
//' Returns a bootstrap (row) by statistic (col) matrix. If \code{return_auc = TRUE},
//' the output has a single column, the AUC. If \code{return_auc = FALSE}, the output
//' has a single column for each evaluation time.
//'  
//' @param boot Bootstrap replicates.
//' @param eval_times Evaluation times.
//' @param idx Unique subject index. 
//' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
//' @param time Observation time.
//' @param value Observation value.
//' @param replace_na Replace NaN with zero?
//' @param return_auc Return the AUC?
//' @param trunc_time Truncation time? Optional. If omitted, defaults
//' to the maximum evaluation time.
//' @return Numeric matrix.
// [[Rcpp::export]]

SEXP BootstrapSamplesR(
  const int boot,
  arma::colvec eval_times,
  const arma::colvec idx,
  const arma::colvec status,
  const arma::colvec time,
  const arma::colvec value,
  const bool replace_na=false,
  const bool return_auc=false,
  const Rcpp::Nullable<double> trunc_time=R_NilValue
){
  
  // Truncation time.
  double tau;
  if(trunc_time.isNotNull()) {
    tau = Rcpp::as<double>(trunc_time);
    eval_times = Truncate(eval_times, tau);
  } else {
    tau = eval_times.max();
  }

  // Structure to store bootstrap samples.
  arma::mat samples;
  if(return_auc) {
    samples = arma::zeros(boot);
  } else {
    const int n_times = eval_times.size();
    samples = arma::zeros(n_times, boot);
  }

  // Loop over bootstrap samples.
  for(int b=0; b<boot; b++) {

    arma::mat boot_data = DrawBootstrapCpp(idx, status, time, value);

    arma::colvec boot_idx = boot_data.col(0);
    arma::colvec boot_status = boot_data.col(1);
    arma::colvec boot_time = boot_data.col(2);
    arma::colvec boot_value = boot_data.col(3);

    arma::mat est = EstimatorCpp(
      eval_times,
      boot_idx,
      boot_status,
      boot_time,
      tau,
      boot_value,
      replace_na, 
      return_auc
    );

    if(return_auc) {
      samples(b, 0) = est(0, 0);
    } else {
      samples.col(b) = est.col(5);
    }
  }

  return Rcpp::wrap(arma::trans(samples));
}


// ----------------------------------------------------------------------------
// Influence function.
// ----------------------------------------------------------------------------


//' Calculate Mu R
//'
//' Evaluate \eqn{\mu(t; \tau) = \int_{t}^{\tau}{S(u)d(u)/y(u)}du}.
//' 
//' @param d Value of d(t) at each time point.
//' @param surv Value of S(t) at each time point.
//' @param unique_times Unique values of time t.
//' @param y Value of y(t) at each time point.
//' @return Numeric vector of \eqn{\mu(t; tau)}.
// [[Rcpp::export]]
SEXP CalcMuR(
  const arma::colvec d,
  const arma::colvec surv,
  const arma::colvec unique_times,
  const arma::colvec y
) {
  
  // Calculate dt.
  const int n_time = unique_times.size();
  arma::colvec delta_t = arma::zeros(n_time);
  delta_t.subvec(1, n_time - 1) = arma::diff(unique_times);

  // Calculate S(t)d(t) / y(t).
  arma::colvec integrand = surv % (d / y);
  
  // Integrate.
  arma::colvec out = arma::reverse(
    arma::cumsum(arma::reverse(delta_t % integrand))
  );
  return Rcpp::wrap(out);
}


// Calculate Mu Cpp
//
// Evaluate \eqn{\mu(t; \tau) = \int_{t}^{\tau}{S(u)d(u)/y(u)}du}.
// 
// @param d Value of d(t) at each time point.
// @param surv Value of S(t) at each time point.
// @param unique_times Unique values of time t.
// @param y Value of y(t) at each time point.
// @return Numeric vector of \eqn{\mu(t; tau)}.
arma::colvec CalcMuCpp(
  const arma::colvec d,
  const arma::colvec surv,
  const arma::colvec unique_times,
  const arma::colvec y
) {
  
  // Calculate dt.
  const int n_time = unique_times.size();
  arma::colvec delta_t = arma::zeros(n_time);
  delta_t.subvec(1, n_time - 1) = arma::diff(unique_times);
  
  // Calculate S(t)d(t) / y(t).
  arma::colvec integrand = surv % (d / y);
  
  // Integrate.
  arma::colvec out = arma::reverse(
    arma::cumsum(arma::reverse(delta_t % integrand))
  );
  return out;
}


// ----------------------------------------------------------------------------


//' Calculate Martingale
//'
//' Calculate \eqn{dM_{i}(t) = dN_{i}(t) - Y_{i}(t)d\Lambda(t)}.
//' 
//' @param haz Value of the hazard at each unique time.
//' @param idx Subject index.
//' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
//' @param time Subject observation times.
//' @param unique_times Unique times at which to obtain the martingale.
//' @return Matrix with subjects as rows and unique times as columns.
// [[Rcpp::export]]
SEXP CalcMartingaleR(
  const arma::colvec haz,
  const arma::colvec idx,
  const arma::colvec status,
  const arma::colvec time,
  const arma::colvec unique_times
) {
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Unique times.
  const int n_times = unique_times.size();
  
  // Create a subject by evaluation times matrix, where
  // dM[i, t] is the martingale increment for subject i at time t.
  arma::mat dm = arma::zeros(n, n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Time and status for the focus subject.
    arma::colvec subj_times = time.elem(arma::find(idx == unique_idx(i)));
    arma::colvec subj_status = status.elem(arma::find(idx == unique_idx(i)));
    
    // Subject is at risk until time > the subject's last time.
    double subj_last_time = subj_times.max();
    
    // Subject's final status.
    arma::uvec key_final_status = arma::find(
      subj_times == subj_last_time, 1, "last");
    double final_status = arma::as_scalar(
      subj_status.elem(key_final_status)
    );
    
    // Loop over times.
    for(int j=0; j<n_times; j++) {
      
      double current_time = unique_times(j);
      double current_haz = haz(j);
      
      // Add dN_{i}(t).
      if(current_time == subj_last_time) {
        if(final_status != 0) {
          dm(i, j) += 1;
        }
      }
      
      // Add -Y_{i}(t)d\Lambda(t).
      if(current_time <= subj_last_time) {
        dm(i, j) -= current_haz;
      } else {
        break;
      }
    }
  }
  return Rcpp::wrap(dm);
}


// Calculate Martingale
//
// Calculate \eqn{dM_{i}(t) = dN_{i}(t) - Y_{i}(t)d\Lambda(t)}.
// 
// @param haz Value of the hazard at each unique time.
// @param idx Subject index.
// @param status Subject status.
// @param time Subject observation times.
// @param unique_times Unique times at which to obtain the martingale.
// @return Matrix with subjects as rows and unique times as columns.
// [[Rcpp::export]]
arma::mat CalcMartingaleCpp(
    const arma::colvec haz,
    const arma::colvec idx,
    const arma::colvec status,
    const arma::colvec time,
    const arma::colvec unique_times
) {
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Unique times.
  const int n_times = unique_times.size();
  
  // Create a subject by evaluation times matrix, where
  // dM[i, t] is the martingale increment for subject i at time t.
  arma::mat dm = arma::zeros(n, n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Time and status for the focus subject.
    arma::colvec subj_times = time.elem(arma::find(idx == unique_idx(i)));
    arma::colvec subj_status = status.elem(arma::find(idx == unique_idx(i)));
    
    // Subject is at risk until time > the subject's las time.
    double subj_last_time = subj_times.max();
    
    // Subject's final status.
    arma::uvec key_final_status = arma::find(
      subj_times == subj_last_time, 1, "last");
    double final_status = arma::as_scalar(
      subj_status.elem(key_final_status)
    );
    
    // Loop over times.
    for(int j=0; j<n_times; j++) {
      
      double current_time = unique_times(j);
      double current_haz = haz(j);
      
      // Add dN_{i}(t).
      if(current_time == subj_last_time) {
        if(final_status != 0) {
          dm(i, j) += 1;
        }
      }
      
      // Add -Y_{i}(t)d\Lambda(t).
      if(current_time <= subj_last_time) {
        dm(i, j) -= current_haz;
      } else {
        break;  
      }
    }
  }
  return dm;
}


// ----------------------------------------------------------------------------


// Calculate I1
//
// Calculate \eqn{I_{1,i} = \int_{0}^{\tau} -\mu(t; \tau) dM_{i}(t) / y(t)}.
// 
// @param dm Matrix of dM_{i}(t).
// @param mu Vector of mu(t; tau).
// @param y Vector of y(t).
// @return Vector with I1 for each subject.
arma::colvec CalcI1Cpp(
    const arma::mat dm,
    const arma::colvec mu,
    const arma::colvec y
) {
  int n = dm.n_rows;
  arma::colvec out = arma::zeros(n);
  for(int i=0; i<n; i++) {
    arma::colvec dmi = arma::trans(dm.row(i));
    out(i) = -1 * arma::sum(mu % dmi / y);
  }
  return out;
}


// Calculate I2
//
// Calculate \eqn{I_{2,i} = \int_{0}^{\tau} S(t)\{D_{i}(t) - d(t)\}/y(t) dt}.
// 
// @param d Vector of d(t).
// @param surv Vector of S(t).
// @param unique_times Vector of unique times t.
// @param value_mat Matrix of D_{i}(t).
// @param y Vector of y(t).
// @return Vector with I2 for each subject.
arma::colvec CalcI2Cpp(
    const arma::mat d,
    const arma::colvec surv,
    const arma::colvec unique_times,
    const arma::mat value_mat,
    const arma::colvec y
) {
  int n = value_mat.n_rows;
  arma::colvec out = arma::zeros(n);
  
  // Calculate dt.
  const int n_time = unique_times.size();
  arma::colvec delta_t = arma::zeros(n_time);
  delta_t.subvec(0, n_time - 2) = arma::diff(unique_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    arma::colvec di = arma::trans(value_mat.row(i));
    
    // Integrate.
    out(i) = arma::sum(surv % (di - d) % delta_t / y);
  }
  return out;
}


// Calculate I3
//
// Calculate \eqn{I_{3,i} = \int_{0}^{\tau} -S(t)d(t)\{Y_{i}(t) - y(t)\} / y^2(t) dt}.
// 
// @param d Vector of d(t).
// @param risk_mat Matrix of Y_{i}(t).
// @param surv Vector of S(t).
// @param unique_times Vector of unique times t.
// @param y Vector of y(t).
// @return Vector with I3 for each subject.
arma::colvec CalcI3Cpp(
    const arma::mat d,
    const arma::mat risk_mat,
    const arma::colvec surv,
    const arma::colvec unique_times,
    const arma::colvec y
) {
  int n = risk_mat.n_rows;
  arma::colvec out = arma::zeros(n);
  
  // Calculate dt.
  const int n_time = unique_times.size();
  arma::colvec delta_t = arma::zeros(n_time);
  delta_t.subvec(0, n_time - 2) = arma::diff(unique_times);
  
  const arma::colvec y2 = arma::pow(y, 2);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    arma::colvec yi = arma::trans(risk_mat.row(i));
    
    // Integrate.
    out(i) = -1 * arma::sum(surv % d % (yi - y) % delta_t / y2);
  }
  return out;
}


// ----------------------------------------------------------------------------


//' Influence Function R
//' 
//' Influence function contributions for the AUC. Includes the three component
//' integrals (`i1`, `i2`, `i3`) and the overall influence `psi`.
//'  
//' @param idx Unique subject index. 
//' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
//' @param time Observation time.
//' @param trunc_time Truncation time? Optional. If omitted, defaults
//' to the maximum evaluation time.
//' @param value Observation value.
//' @return Data.frame.
// [[Rcpp::export]]

SEXP InfluenceR(
    const arma::colvec idx,
    const arma::colvec status,
    const arma::colvec time,
    const double trunc_time,
    const arma::colvec value
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  // const int n = unique_idx.size();
  
  // Unique times.
  arma::colvec unique_times = arma::unique(time);
  unique_times = Truncate(unique_times, trunc_time);
  // const int n_times = unique_times.size();
  
  // Return columns: {0: time, 1: y, 2: haz, 3: surv, 4: d, 5: exp}.
  arma::mat est = EstimatorCpp(
    unique_times, idx, status, time, trunc_time, value);
  // Rcpp::Rcout << est << std::endl; 
  
  const arma::colvec y = est.col(1);
  const arma::colvec haz = est.col(2);
  const arma::colvec surv = est.col(3);
  const arma::colvec d = est.col(4);
  
  // Calculate mu.
  const arma::colvec mu = CalcMuCpp(d, surv, unique_times, y);
  // Rcpp::Rcout << mu << std::endl; 
  
  // Calculate martingales.
  const arma::mat dm = CalcMartingaleCpp(haz, idx, status, time, unique_times);
  // Rcpp::Rcout << dm << std::endl; 
  
  // Calculate I1.
  const arma::colvec i1 = CalcI1Cpp(dm, mu, y);
  // Rcpp::Rcout << i1 << std::endl; 
  
  // Calculate value matrix.
  const arma::mat value_mat = ValueMatrixCpp(unique_times, idx, time, value);
  
  // Calculate I2.
  const arma::colvec i2 = CalcI2Cpp(d, surv, unique_times, value_mat, y);
  // Rcpp::Rcout << i2 << std::endl; 
  
  // Calculate at risk matrix.
  const arma::mat risk_mat = AtRiskMatrixCpp(unique_times, idx, time);
  
  // Calculate I3.
  const arma::colvec i3 = CalcI3Cpp(d, risk_mat, surv, unique_times, y);
  // Rcpp::Rcout << i3 << std::endl; 
  
  // Overall influence funciton.
  const arma::colvec psi = i1 + i2 + i3;
  
  // Output.
  return Rcpp::DataFrame::create(
    Rcpp::Named("idx")=unique_idx,
    Rcpp::Named("i1")=i1,
    Rcpp::Named("i2")=i2,
    Rcpp::Named("i3")=i3,
    Rcpp::Named("psi")=psi
  );
}


// Influence Function R
// 
// Influence function contributions for the AUC. Returns the overall influence
// only.
// 
// @param idx Unique subject index. 
// @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event. 
// @param time Observation time.
// @param trunc_time Truncation time? Optional. If omitted, defaults
// to the maximum evaluation time.
// @param value Observation value.
// @return Numeric vector.

arma::colvec InfluenceCpp(
    const arma::colvec idx,
    const arma::colvec status,
    const arma::colvec time,
    const double trunc_time,
    const arma::colvec value
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  // const int n = unique_idx.size();
  
  // Unique times.
  arma::colvec unique_times = arma::unique(time);
  unique_times = Truncate(unique_times, trunc_time);
  // const int n_times = unique_times.size();
  
  // Return columns: {0: time, 1: y, 2: haz, 3: surv, 4: d, 5: exp}.
  arma::mat est = EstimatorCpp(
    unique_times, idx, status, time, trunc_time, value);
  // Rcpp::Rcout << est << std::endl; 
  
  const arma::colvec y = est.col(1);
  const arma::colvec haz = est.col(2);
  const arma::colvec surv = est.col(3);
  const arma::colvec d = est.col(4);
  
  // Calculate mu.
  const arma::colvec mu = CalcMuCpp(d, surv, unique_times, y);
  // Rcpp::Rcout << mu << std::endl; 
  
  // Calculate martingales.
  const arma::mat dm = CalcMartingaleCpp(haz, idx, status, time, unique_times);
  // Rcpp::Rcout << dm << std::endl; 
  
  // Calculate I1.
  const arma::colvec i1 = CalcI1Cpp(dm, mu, y);
  // Rcpp::Rcout << i1 << std::endl; 
  
  // Calculate value matrix.
  const arma::mat value_mat = ValueMatrixCpp(unique_times, idx, time, value);
  
  // Calculate I2.
  const arma::colvec i2 = CalcI2Cpp(d, surv, unique_times, value_mat, y);
  // Rcpp::Rcout << i2 << std::endl; 
  
  // Calculate at risk matrix.
  const arma::mat risk_mat = AtRiskMatrixCpp(unique_times, idx, time);
  
  // Calculate I3.
  const arma::colvec i3 = CalcI3Cpp(d, risk_mat, surv, unique_times, y);
  // Rcpp::Rcout << i3 << std::endl; 
  
  // Overall influence funciton.
  const arma::colvec psi = i1 + i2 + i3;
  
  // Output.
  return psi;
}


// ----------------------------------------------------------------------------

//' Perturbation R
//'  
//' Generates realizations of \eqn{\frac{1}{n}\sum_{i=1}^{n}\psi_{i}w_{i}},
//' where \eqn{\psi_{i}} is the influence function for the ith subject and the
//' \eqn{w_{i}} are IID random weights.
//' 
//' @section Notes:
//' The random seed should be set in R prior to calling this function.
//'  
//' @param idx Unique subject index. 
//' @param perturbations Number of perturbations
//' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
//' @param time Observation time.
//' @param trunc_time Truncation time? Optional. If omitted, defaults
//' to the maximum evaluation time.
//' @param value Observation value.
//' @return Numeric vector.
// [[Rcpp::export]]

SEXP PerturbationR(
    const arma::colvec idx,
    const int perturbations,
    const arma::colvec status,
    const arma::colvec time,
    const double trunc_time,
    const arma::colvec value
){

  arma::colvec out = arma::zeros(perturbations);
  arma::colvec psi = InfluenceCpp(idx, status, time, trunc_time, value);
  arma::colvec weights(psi.size());
  
  for(int i=0; i<perturbations; i++) {
    weights = arma::randn(psi.size());
    out(i) = arma::mean(psi % weights);
  }
  
  return Rcpp::wrap(out);
}


// ----------------------------------------------------------------------------


//' Interpolation R
//'  
//' Linearly interpolations between each subject's measurements.
//' The input data should contain no missing values. 
//'  
//' @param idx Unique subject index. 
//' @param status Status, coded as 0 for censoring, 1 for event, 2 for terminal event.
//' @param time Observation time.
//' @param value Observation value.
//' @param n_points Number of interpolation points.
//' @return Data.frame.
//' @export
// [[Rcpp::export]]

SEXP InterpolateR(
    const arma::colvec idx,
    const arma::colvec status,
    const arma::colvec time,
    const arma::colvec value,
    const int n_points = 100
){
  
  // Subjects.
  const arma::colvec unique_idx = arma::unique(idx);
  const int n = unique_idx.size();
  
  // Output structure.
  arma::colvec out_idx = arma::zeros(n * n_points);
  arma::colvec out_status = arma::ones(n * n_points);
  arma::colvec out_time = arma::zeros(n * n_points);
  arma::colvec out_value = arma::zeros(n * n_points);

  for(int i=0; i<n; i++) {
    
    // Current subject.
    arma::colvec subj_status = status.elem(arma::find(idx == unique_idx(i)));
    double final_status = subj_status(subj_status.size() - 1);
    arma::colvec subj_time = time.elem(arma::find(idx == unique_idx(i)));
    arma::colvec subj_value = value.elem(arma::find(idx == unique_idx(i)));
    
    // Interpolation points.
    arma::colvec grid = arma::linspace<arma::colvec>(0, subj_time.max(), n_points);
    arma::colvec yhat;
    arma::interp1(subj_time, subj_value, grid, yhat, "linear");
    
    // Store.
    int first_idx = i * n_points;
    int last_idx = (i + 1) * n_points - 1;
    out_idx(arma::span(first_idx, last_idx)) = unique_idx(i) * arma::ones(n_points);
    out_status(last_idx) = final_status;
    out_time(arma::span(first_idx, last_idx)) = grid;
    out_value(arma::span(first_idx, last_idx)) = yhat;
  }
  
  // Output.
  return Rcpp::DataFrame::create(
    Rcpp::Named("idx")=out_idx,
    Rcpp::Named("status")=out_status,
    Rcpp::Named("time")=out_time,
    Rcpp::Named("value")=out_value
  );
}