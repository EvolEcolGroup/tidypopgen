// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fbm256_prod_and_rowSumsSq
List fbm256_prod_and_rowSumsSq(Environment BM, const IntegerVector& ind_row, const IntegerVector& ind_col, const NumericVector& center, const NumericVector& scale, const NumericMatrix& V);
RcppExport SEXP _tidypopgen_fbm256_prod_and_rowSumsSq(SEXP BMSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(fbm256_prod_and_rowSumsSq(BM, ind_row, ind_col, center, scale, V));
    return rcpp_result_gen;
END_RCPP
}
// gt_grouped_alt_freq_diploid
ListOf<NumericMatrix> gt_grouped_alt_freq_diploid(Environment BM, const IntegerVector& rowInd, const IntegerVector& colInd, const IntegerVector& groupIds, int ngroups, int ncores);
RcppExport SEXP _tidypopgen_gt_grouped_alt_freq_diploid(SEXP BMSEXP, SEXP rowIndSEXP, SEXP colIndSEXP, SEXP groupIdsSEXP, SEXP ngroupsSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type groupIds(groupIdsSEXP);
    Rcpp::traits::input_parameter< int >::type ngroups(ngroupsSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(gt_grouped_alt_freq_diploid(BM, rowInd, colInd, groupIds, ngroups, ncores));
    return rcpp_result_gen;
END_RCPP
}
// gt_grouped_alt_freq_pseudohap
ListOf<NumericMatrix> gt_grouped_alt_freq_pseudohap(Environment BM, const IntegerVector& rowInd, const IntegerVector& colInd, const IntegerVector& groupIds, int ngroups, const IntegerVector& ploidy, int ncores);
RcppExport SEXP _tidypopgen_gt_grouped_alt_freq_pseudohap(SEXP BMSEXP, SEXP rowIndSEXP, SEXP colIndSEXP, SEXP groupIdsSEXP, SEXP ngroupsSEXP, SEXP ploidySEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type groupIds(groupIdsSEXP);
    Rcpp::traits::input_parameter< int >::type ngroups(ngroupsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(gt_grouped_alt_freq_pseudohap(BM, rowInd, colInd, groupIds, ngroups, ploidy, ncores));
    return rcpp_result_gen;
END_RCPP
}
// gt_grouped_missingness
NumericMatrix gt_grouped_missingness(Environment BM, const IntegerVector& rowInd, const IntegerVector& colInd, const IntegerVector& groupIds, int ngroups, int ncores);
RcppExport SEXP _tidypopgen_gt_grouped_missingness(SEXP BMSEXP, SEXP rowIndSEXP, SEXP colIndSEXP, SEXP groupIdsSEXP, SEXP ngroupsSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type groupIds(groupIdsSEXP);
    Rcpp::traits::input_parameter< int >::type ngroups(ngroupsSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(gt_grouped_missingness(BM, rowInd, colInd, groupIds, ngroups, ncores));
    return rcpp_result_gen;
END_RCPP
}
// gt_grouped_summaries
ListOf<NumericMatrix> gt_grouped_summaries(Environment BM, const IntegerVector& rowInd, const IntegerVector& colInd, const IntegerVector& groupIds, int ngroups, int ncores);
RcppExport SEXP _tidypopgen_gt_grouped_summaries(SEXP BMSEXP, SEXP rowIndSEXP, SEXP colIndSEXP, SEXP groupIdsSEXP, SEXP ngroupsSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type groupIds(groupIdsSEXP);
    Rcpp::traits::input_parameter< int >::type ngroups(ngroupsSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(gt_grouped_summaries(BM, rowInd, colInd, groupIds, ngroups, ncores));
    return rcpp_result_gen;
END_RCPP
}
// SNPHWE2
double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp);
RcppExport SEXP _tidypopgen_SNPHWE2(SEXP obs_hetsSEXP, SEXP obs_hom1SEXP, SEXP obs_hom2SEXP, SEXP midpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int32_t >::type obs_hets(obs_hetsSEXP);
    Rcpp::traits::input_parameter< int32_t >::type obs_hom1(obs_hom1SEXP);
    Rcpp::traits::input_parameter< int32_t >::type obs_hom2(obs_hom2SEXP);
    Rcpp::traits::input_parameter< uint32_t >::type midp(midpSEXP);
    rcpp_result_gen = Rcpp::wrap(SNPHWE2(obs_hets, obs_hom1, obs_hom2, midp));
    return rcpp_result_gen;
END_RCPP
}
// SNPHWE_t
int32_t SNPHWE_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);
RcppExport SEXP _tidypopgen_SNPHWE_t(SEXP obs_hetsSEXP, SEXP obs_hom1SEXP, SEXP obs_hom2SEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int32_t >::type obs_hets(obs_hetsSEXP);
    Rcpp::traits::input_parameter< int32_t >::type obs_hom1(obs_hom1SEXP);
    Rcpp::traits::input_parameter< int32_t >::type obs_hom2(obs_hom2SEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(SNPHWE_t(obs_hets, obs_hom1, obs_hom2, thresh));
    return rcpp_result_gen;
END_RCPP
}
// SNPHWE_midp_t
int32_t SNPHWE_midp_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);
RcppExport SEXP _tidypopgen_SNPHWE_midp_t(SEXP obs_hetsSEXP, SEXP obs_hom1SEXP, SEXP obs_hom2SEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int32_t >::type obs_hets(obs_hetsSEXP);
    Rcpp::traits::input_parameter< int32_t >::type obs_hom1(obs_hom1SEXP);
    Rcpp::traits::input_parameter< int32_t >::type obs_hom2(obs_hom2SEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(SNPHWE_midp_t(obs_hets, obs_hom1, obs_hom2, thresh));
    return rcpp_result_gen;
END_RCPP
}
// increment_as_counts
void increment_as_counts(Environment k, Environment k2, arma::mat& na_mat, arma::mat& dos_mat, Environment BM, const IntegerVector& rowInd, const IntegerVector& colInd);
RcppExport SEXP _tidypopgen_increment_as_counts(SEXP kSEXP, SEXP k2SEXP, SEXP na_matSEXP, SEXP dos_matSEXP, SEXP BMSEXP, SEXP rowIndSEXP, SEXP colIndSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type k(kSEXP);
    Rcpp::traits::input_parameter< Environment >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type na_mat(na_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type dos_mat(dos_matSEXP);
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    increment_as_counts(k, k2, na_mat, dos_mat, BM, rowInd, colInd);
    return R_NilValue;
END_RCPP
}
// increment_ibs_counts
void increment_ibs_counts(Environment k, Environment k2, arma::mat& genotype0, arma::mat& genotype1, arma::mat& genotype2, Environment BM, const IntegerVector& rowInd, const IntegerVector& colInd);
RcppExport SEXP _tidypopgen_increment_ibs_counts(SEXP kSEXP, SEXP k2SEXP, SEXP genotype0SEXP, SEXP genotype1SEXP, SEXP genotype2SEXP, SEXP BMSEXP, SEXP rowIndSEXP, SEXP colIndSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type k(kSEXP);
    Rcpp::traits::input_parameter< Environment >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type genotype0(genotype0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type genotype1(genotype1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type genotype2(genotype2SEXP);
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    increment_ibs_counts(k, k2, genotype0, genotype1, genotype2, BM, rowInd, colInd);
    return R_NilValue;
END_RCPP
}
// increment_king_numerator
void increment_king_numerator(Environment k, Environment n_Aa_i, arma::mat& genotype0, arma::mat& genotype1, arma::mat& genotype2, arma::mat& genotype_valid, Environment BM, const IntegerVector& rowInd, const IntegerVector& colInd);
RcppExport SEXP _tidypopgen_increment_king_numerator(SEXP kSEXP, SEXP n_Aa_iSEXP, SEXP genotype0SEXP, SEXP genotype1SEXP, SEXP genotype2SEXP, SEXP genotype_validSEXP, SEXP BMSEXP, SEXP rowIndSEXP, SEXP colIndSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type k(kSEXP);
    Rcpp::traits::input_parameter< Environment >::type n_Aa_i(n_Aa_iSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type genotype0(genotype0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type genotype1(genotype1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type genotype2(genotype2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type genotype_valid(genotype_validSEXP);
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    increment_king_numerator(k, n_Aa_i, genotype0, genotype1, genotype2, genotype_valid, BM, rowInd, colInd);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tidypopgen_fbm256_prod_and_rowSumsSq", (DL_FUNC) &_tidypopgen_fbm256_prod_and_rowSumsSq, 6},
    {"_tidypopgen_gt_grouped_alt_freq_diploid", (DL_FUNC) &_tidypopgen_gt_grouped_alt_freq_diploid, 6},
    {"_tidypopgen_gt_grouped_alt_freq_pseudohap", (DL_FUNC) &_tidypopgen_gt_grouped_alt_freq_pseudohap, 7},
    {"_tidypopgen_gt_grouped_missingness", (DL_FUNC) &_tidypopgen_gt_grouped_missingness, 6},
    {"_tidypopgen_gt_grouped_summaries", (DL_FUNC) &_tidypopgen_gt_grouped_summaries, 6},
    {"_tidypopgen_SNPHWE2", (DL_FUNC) &_tidypopgen_SNPHWE2, 4},
    {"_tidypopgen_SNPHWE_t", (DL_FUNC) &_tidypopgen_SNPHWE_t, 4},
    {"_tidypopgen_SNPHWE_midp_t", (DL_FUNC) &_tidypopgen_SNPHWE_midp_t, 4},
    {"_tidypopgen_increment_as_counts", (DL_FUNC) &_tidypopgen_increment_as_counts, 7},
    {"_tidypopgen_increment_ibs_counts", (DL_FUNC) &_tidypopgen_increment_ibs_counts, 8},
    {"_tidypopgen_increment_king_numerator", (DL_FUNC) &_tidypopgen_increment_king_numerator, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_tidypopgen(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
