/******************************************************************************/

#include "fbm_byte_acc.h"

/******************************************************************************/

// Grouped (per-population) version of gt_pi_polyploid(); see that function
// for details of the estimator. When the dataset is marked pseudohaploid,
// individuals with ploidy one retain 0/2 storage while diploid individuals
// retain their usual dosage. The dataset-level flag distinguishes this from
// ordinary 0/1 haploid calls. The number of valid allele copies (`n`) is
// only returned if return_n = TRUE.

// [[Rcpp::export]]
ListOf<NumericMatrix> gt_grouped_pi_polyploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   int ngroups,
                                   const NumericVector& ploidy,
                                   bool is_pseudohaploid,
                                   int ncores,
                                   int block_size = 100000,
                                   bool return_n = false) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }
  std::vector<size_t> group_ids = check_group_ids(groupIds, n, ngroups);

  // pseudohaploid individuals contribute half their (0/2) dosage, and every
  // individual contributes as many valid alleles as its ploidy
  std::vector<double> weight(n), copies(n);
  for (size_t i = 0; i < n; i++) {
    weight[i] = (is_pseudohaploid && ploidy[i] == 1) ? 0.5 : 1;
    copies[i] = ploidy[i];
  }

  NumericMatrix pi(m, ngroups);
  // only allocate the valid allele counts if they are needed
  NumericMatrix valid_alleles(return_n ? m : 0, return_n ? ngroups : 0);
  double* p_pi = pi.begin();
  double* p_valid = valid_alleles.begin();

  alt_valid_counts(acc, group_ids, ngroups, weight, copies, ncores, block_size,
                   [&](size_t j, size_t g, double alt,
                       double this_valid_alleles) {
                     // guard against 0 or 1 valid allele copies in this group
                     // at this locus, where the denominator would be 0
                     p_pi[j + g * m] = (this_valid_alleles > 1) ?
                       (alt * (this_valid_alleles - alt) /
                         (this_valid_alleles * (this_valid_alleles - 1) / 2)) :
                       NA_REAL;
                     if (return_n) p_valid[j + g * m] = this_valid_alleles;
                   });

  if (!return_n) {
    return List::create(_["pi"] = pi);
  }
  return List::create(_["pi"]  = pi,
                      _["n"] = valid_alleles);
}

/******************************************************************************/
