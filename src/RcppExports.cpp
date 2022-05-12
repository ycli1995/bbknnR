// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// UpdateNumMat
NumericMatrix UpdateNumMat(NumericMatrix& raw_mat, NumericMatrix& new_mat, IntegerVector& row_idx, IntegerVector& col_idx);
RcppExport SEXP _bbknnR_UpdateNumMat(SEXP raw_matSEXP, SEXP new_matSEXP, SEXP row_idxSEXP, SEXP col_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type raw_mat(raw_matSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type new_mat(new_matSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type row_idx(row_idxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type col_idx(col_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateNumMat(raw_mat, new_mat, row_idx, col_idx));
    return rcpp_result_gen;
END_RCPP
}
// UpdateIntMat
IntegerMatrix UpdateIntMat(IntegerMatrix& raw_mat, IntegerMatrix& new_mat, IntegerVector& row_idx, IntegerVector& col_idx);
RcppExport SEXP _bbknnR_UpdateIntMat(SEXP raw_matSEXP, SEXP new_matSEXP, SEXP row_idxSEXP, SEXP col_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix& >::type raw_mat(raw_matSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type new_mat(new_matSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type row_idx(row_idxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type col_idx(col_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateIntMat(raw_mat, new_mat, row_idx, col_idx));
    return rcpp_result_gen;
END_RCPP
}
// GetRawIndex
IntegerMatrix GetRawIndex(IntegerMatrix& idx_mat, IntegerVector& raw_idx);
RcppExport SEXP _bbknnR_GetRawIndex(SEXP idx_matSEXP, SEXP raw_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix& >::type idx_mat(idx_matSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type raw_idx(raw_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(GetRawIndex(idx_mat, raw_idx));
    return rcpp_result_gen;
END_RCPP
}
// GetSparseDist
List GetSparseDist(IntegerMatrix& knn_index, NumericMatrix& knn_dist, int n_obs, int n_neighbors);
RcppExport SEXP _bbknnR_GetSparseDist(SEXP knn_indexSEXP, SEXP knn_distSEXP, SEXP n_obsSEXP, SEXP n_neighborsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix& >::type knn_index(knn_indexSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type knn_dist(knn_distSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< int >::type n_neighbors(n_neighborsSEXP);
    rcpp_result_gen = Rcpp::wrap(GetSparseDist(knn_index, knn_dist, n_obs, n_neighbors));
    return rcpp_result_gen;
END_RCPP
}
// Trimming
void Trimming(NumericVector& x, IntegerVector& row_idx, IntegerVector& p, int trim);
RcppExport SEXP _bbknnR_Trimming(SEXP xSEXP, SEXP row_idxSEXP, SEXP pSEXP, SEXP trimSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type row_idx(row_idxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type trim(trimSEXP);
    Trimming(x, row_idx, p, trim);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bbknnR_UpdateNumMat", (DL_FUNC) &_bbknnR_UpdateNumMat, 4},
    {"_bbknnR_UpdateIntMat", (DL_FUNC) &_bbknnR_UpdateIntMat, 4},
    {"_bbknnR_GetRawIndex", (DL_FUNC) &_bbknnR_GetRawIndex, 2},
    {"_bbknnR_GetSparseDist", (DL_FUNC) &_bbknnR_GetSparseDist, 4},
    {"_bbknnR_Trimming", (DL_FUNC) &_bbknnR_Trimming, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_bbknnR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}