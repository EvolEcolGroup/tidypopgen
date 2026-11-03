/******************************************************************************/

#include "fbm_byte_acc.h"

/******************************************************************************/

// lookup table flagging the bytes that decode to a missing genotype
static std::vector<uint64_t> missing_lut(const FBMByteAcc& acc) {
  std::vector<uint64_t> lut(256);
  for (int b = 0; b < 256; b++) {
    lut[b] = (acc.code(b) > -1) ? 0 : 1;
  }
  return lut;
}

/******************************************************************************/

// Count the missing genotypes at each locus for each group or, if as_counts =
// FALSE, compute the proportion of missing genotypes (dividing by the number
// of individuals in the group; NaN for empty groups). Loci are processed in
// blocks of block_size, checking for user interrupts between blocks.

// [[Rcpp::export]]
NumericMatrix grouped_missingness_cpp(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   int ngroups,
                                   int ncores,
                                   bool as_counts = true,
                                   int block_size = 100000) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  std::vector<size_t> group_ids = check_group_ids(groupIds, n, ngroups);
  // number of individuals in each group (to compute proportions)
  std::vector<double> group_sizes(ngroups, 0);
  for (size_t i = 0; i < n; i++) group_sizes[group_ids[i]] += 1;

  std::vector<uint64_t> lut = missing_lut(acc);
  // all individuals use the same lookup table
  GroupedRows grouped(acc.rows(), group_ids, std::vector<size_t>(n, 0),
                      max_run_length<1>(lut));

  NumericMatrix na_counts(m, ngroups);
  double* out = na_counts.begin();

  lut_sums<1>(acc, grouped, lut, ngroups, ncores, block_size,
              [&](size_t j, const uint64_t* tot) {
                for (int g = 0; g < ngroups; g++) {
                  double n_na = static_cast<double>(tot[g]);
                  out[j + g * m] = as_counts ? n_na : n_na / group_sizes[g];
                }
              });
  return na_counts;
}

/******************************************************************************/

// Count the missing genotypes at each locus (all individuals together)

// [[Rcpp::export]]
IntegerVector missingness_cpp(Environment BM,
                              const IntegerVector& rowInd,
                              const IntegerVector& colInd,
                              int ncores,
                              int block_size = 100000) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  std::vector<uint64_t> lut = missing_lut(acc);
  // all individuals belong to one group and use the same lookup table
  GroupedRows grouped(acc.rows(), std::vector<size_t>(n, 0),
                      std::vector<size_t>(n, 0), max_run_length<1>(lut));

  IntegerVector n_na(m);
  int* out = n_na.begin();

  lut_sums<1>(acc, grouped, lut, 1, ncores, block_size,
              [&](size_t j, const uint64_t* tot) {
                out[j] = static_cast<int>(tot[0]);
              });
  return n_na;
}

/******************************************************************************/
