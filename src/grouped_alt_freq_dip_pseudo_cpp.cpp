/******************************************************************************/

#include "fbm_byte_acc.h"

/******************************************************************************/

// Estimate the alternate allele frequency at each locus for each group, for
// diploid and pseudohaploid data. If freq_only = TRUE (and as_counts = FALSE),
// only the ngroups columns of frequencies are returned, to avoid allocating
// the valid allele counts when they are not needed. Loci are processed in
// blocks of block_size, checking for user interrupts between blocks.

// [[Rcpp::export]]
NumericMatrix grouped_alt_freq_dip_pseudo_cpp(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   int ngroups,
                                   const NumericVector& ploidy,
                                   int ncores,
                                   bool as_counts,
                                   int block_size = 100000,
                                   bool freq_only = false) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }
  std::vector<size_t> group_ids = check_group_ids(groupIds, n, ngroups);

  // multiplier to convert dosages for pseudohaploid data, and number of
  // valid alleles contributed by each individual
  std::vector<double> pseudo_multiplier(n), copies(n);
  for (size_t i = 0; i < n; i++) {
    pseudo_multiplier[i] = 1 / (3 - ploidy[i]);
    copies[i] = ploidy[i];
  }

  // We use a single large matrix, the first ngroup columns are frequencies
  // and the next ngroups columns are valid alleles (unless freq_only). If
  // as_counts = TRUE, then we return the counts instead of the frequencies
  bool with_valid = as_counts || !freq_only;
  NumericMatrix freq_mat(m, with_valid ? ngroups * 2 : ngroups);
  double* out = freq_mat.begin();

  alt_valid_counts(acc, group_ids, ngroups, pseudo_multiplier, copies, ncores,
                   block_size,
                   [&](size_t j, size_t g, double alt, double valid) {
                     // divide alt count by valid count (NaN if no valid
                     // alleles), unless we want the counts
                     out[j + g * m] = as_counts ? alt : alt / valid;
                     if (with_valid) out[j + (ngroups + g) * m] = valid;
                   });

  return freq_mat;
}

/******************************************************************************/
