// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dir_mul_sample_cpp
IntegerMatrix dir_mul_sample_cpp(int& x, int folds, double b);
RcppExport SEXP _countsplit_dir_mul_sample_cpp(SEXP xSEXP, SEXP foldsSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type folds(foldsSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(dir_mul_sample_cpp(x, folds, b));
    return rcpp_result_gen;
END_RCPP
}
// mapply_dir_mul_sample_cpp
IntegerMatrix mapply_dir_mul_sample_cpp(IntegerVector x, int folds, NumericVector overdisps);
RcppExport SEXP _countsplit_mapply_dir_mul_sample_cpp(SEXP xSEXP, SEXP foldsSEXP, SEXP overdispsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type folds(foldsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type overdisps(overdispsSEXP);
    rcpp_result_gen = Rcpp::wrap(mapply_dir_mul_sample_cpp(x, folds, overdisps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_countsplit_dir_mul_sample_cpp", (DL_FUNC) &_countsplit_dir_mul_sample_cpp, 3},
    {"_countsplit_mapply_dir_mul_sample_cpp", (DL_FUNC) &_countsplit_mapply_dir_mul_sample_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_countsplit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}