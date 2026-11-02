/******************************************************************************/

#include "fbm_byte_acc.h"

/******************************************************************************/

// Estimate nucleotide diversity (pi) at each locus for each group, for
// diploid data. The number of valid alleles at each locus for each group
// (`n`) is only returned if return_n = TRUE (it is not needed by loci_pi()
// and pop_tajimas_d(), and is as large as pi). Loci are processed in blocks
// of block_size, checking for user interrupts between blocks.

// [[Rcpp::export]]
ListOf<NumericMatrix> gt_grouped_pi_diploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   int ngroups,
                                   int ncores,
                                   int block_size = 100000,
                                   bool return_n = false) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  std::vector<size_t> group_ids = check_group_ids(groupIds, n, ngroups);
  // every individual contributes its dosage and two valid alleles
  std::vector<double> weight(n, 1), copies(n, 2);

  NumericMatrix pi(m, ngroups);
  // only allocate the valid allele counts if they are needed
  NumericMatrix valid_alleles(return_n ? m : 0, return_n ? ngroups : 0);
  double* p_pi = pi.begin();
  double* p_valid = valid_alleles.begin();

  alt_valid_counts(acc, group_ids, ngroups, weight, copies, ncores, block_size,
                   [&](size_t j, size_t g, double alt,
                       double this_valid_alleles) {
                     // guard against 0 valid individuals in this group at
                     // this locus, where the denominator would be 0
                     p_pi[j + g * m] = (this_valid_alleles > 0) ?
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
