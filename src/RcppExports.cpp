// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ValueMatrixR
SEXP ValueMatrixR(const arma::colvec eval_times, const arma::colvec idx, const arma::colvec time, const arma::colvec value);
RcppExport SEXP _AURMC_ValueMatrixR(SEXP eval_timesSEXP, SEXP idxSEXP, SEXP timeSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type eval_times(eval_timesSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(ValueMatrixR(eval_times, idx, time, value));
    return rcpp_result_gen;
END_RCPP
}
// AtRiskMatrixR
SEXP AtRiskMatrixR(const arma::colvec eval_times, const arma::colvec idx, const arma::colvec time);
RcppExport SEXP _AURMC_AtRiskMatrixR(SEXP eval_timesSEXP, SEXP idxSEXP, SEXP timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type eval_times(eval_timesSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    rcpp_result_gen = Rcpp::wrap(AtRiskMatrixR(eval_times, idx, time));
    return rcpp_result_gen;
END_RCPP
}
// KaplanMeierR
SEXP KaplanMeierR(const arma::colvec eval_times, const arma::colvec idx, const arma::colvec status, const arma::colvec time);
RcppExport SEXP _AURMC_KaplanMeierR(SEXP eval_timesSEXP, SEXP idxSEXP, SEXP statusSEXP, SEXP timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type eval_times(eval_timesSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    rcpp_result_gen = Rcpp::wrap(KaplanMeierR(eval_times, idx, status, time));
    return rcpp_result_gen;
END_RCPP
}
// EstimatorR
SEXP EstimatorR(const arma::colvec idx, const arma::colvec status, const arma::colvec time, const arma::colvec value, const Rcpp::Nullable<Rcpp::NumericVector> eval_times, const bool replace_na, const bool return_auc, const Rcpp::Nullable<double> trunc_time);
RcppExport SEXP _AURMC_EstimatorR(SEXP idxSEXP, SEXP statusSEXP, SEXP timeSEXP, SEXP valueSEXP, SEXP eval_timesSEXP, SEXP replace_naSEXP, SEXP return_aucSEXP, SEXP trunc_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type value(valueSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericVector> >::type eval_times(eval_timesSEXP);
    Rcpp::traits::input_parameter< const bool >::type replace_na(replace_naSEXP);
    Rcpp::traits::input_parameter< const bool >::type return_auc(return_aucSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<double> >::type trunc_time(trunc_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(EstimatorR(idx, status, time, value, eval_times, replace_na, return_auc, trunc_time));
    return rcpp_result_gen;
END_RCPP
}
// DrawBootstrapR
SEXP DrawBootstrapR(const arma::colvec idx, const arma::colvec status, const arma::colvec time, const arma::colvec value);
RcppExport SEXP _AURMC_DrawBootstrapR(SEXP idxSEXP, SEXP statusSEXP, SEXP timeSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(DrawBootstrapR(idx, status, time, value));
    return rcpp_result_gen;
END_RCPP
}
// BootstrapSamplesR
SEXP BootstrapSamplesR(const int boot, arma::colvec eval_times, const arma::colvec idx, const arma::colvec status, const arma::colvec time, const arma::colvec value, const bool replace_na, const bool return_auc, const Rcpp::Nullable<double> trunc_time);
RcppExport SEXP _AURMC_BootstrapSamplesR(SEXP bootSEXP, SEXP eval_timesSEXP, SEXP idxSEXP, SEXP statusSEXP, SEXP timeSEXP, SEXP valueSEXP, SEXP replace_naSEXP, SEXP return_aucSEXP, SEXP trunc_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type boot(bootSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type eval_times(eval_timesSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type value(valueSEXP);
    Rcpp::traits::input_parameter< const bool >::type replace_na(replace_naSEXP);
    Rcpp::traits::input_parameter< const bool >::type return_auc(return_aucSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<double> >::type trunc_time(trunc_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(BootstrapSamplesR(boot, eval_times, idx, status, time, value, replace_na, return_auc, trunc_time));
    return rcpp_result_gen;
END_RCPP
}
// CalcMuR
SEXP CalcMuR(const arma::colvec d, const arma::colvec surv, const arma::colvec unique_times, const arma::colvec y);
RcppExport SEXP _AURMC_CalcMuR(SEXP dSEXP, SEXP survSEXP, SEXP unique_timesSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type surv(survSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type unique_times(unique_timesSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(CalcMuR(d, surv, unique_times, y));
    return rcpp_result_gen;
END_RCPP
}
// CalcMartingaleR
SEXP CalcMartingaleR(const arma::colvec haz, const arma::colvec idx, const arma::colvec status, const arma::colvec time, const arma::colvec unique_times);
RcppExport SEXP _AURMC_CalcMartingaleR(SEXP hazSEXP, SEXP idxSEXP, SEXP statusSEXP, SEXP timeSEXP, SEXP unique_timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type haz(hazSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type unique_times(unique_timesSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcMartingaleR(haz, idx, status, time, unique_times));
    return rcpp_result_gen;
END_RCPP
}
// CalcMartingaleCpp
arma::mat CalcMartingaleCpp(const arma::colvec haz, const arma::colvec idx, const arma::colvec status, const arma::colvec time, const arma::colvec unique_times);
RcppExport SEXP _AURMC_CalcMartingaleCpp(SEXP hazSEXP, SEXP idxSEXP, SEXP statusSEXP, SEXP timeSEXP, SEXP unique_timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type haz(hazSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type unique_times(unique_timesSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcMartingaleCpp(haz, idx, status, time, unique_times));
    return rcpp_result_gen;
END_RCPP
}
// InfluenceR
SEXP InfluenceR(const arma::colvec idx, const arma::colvec status, const arma::colvec time, const double trunc_time, const arma::colvec value);
RcppExport SEXP _AURMC_InfluenceR(SEXP idxSEXP, SEXP statusSEXP, SEXP timeSEXP, SEXP trunc_timeSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double >::type trunc_time(trunc_timeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(InfluenceR(idx, status, time, trunc_time, value));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_AURMC_ValueMatrixR", (DL_FUNC) &_AURMC_ValueMatrixR, 4},
    {"_AURMC_AtRiskMatrixR", (DL_FUNC) &_AURMC_AtRiskMatrixR, 3},
    {"_AURMC_KaplanMeierR", (DL_FUNC) &_AURMC_KaplanMeierR, 4},
    {"_AURMC_EstimatorR", (DL_FUNC) &_AURMC_EstimatorR, 8},
    {"_AURMC_DrawBootstrapR", (DL_FUNC) &_AURMC_DrawBootstrapR, 4},
    {"_AURMC_BootstrapSamplesR", (DL_FUNC) &_AURMC_BootstrapSamplesR, 9},
    {"_AURMC_CalcMuR", (DL_FUNC) &_AURMC_CalcMuR, 4},
    {"_AURMC_CalcMartingaleR", (DL_FUNC) &_AURMC_CalcMartingaleR, 5},
    {"_AURMC_CalcMartingaleCpp", (DL_FUNC) &_AURMC_CalcMartingaleCpp, 5},
    {"_AURMC_InfluenceR", (DL_FUNC) &_AURMC_InfluenceR, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_AURMC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
