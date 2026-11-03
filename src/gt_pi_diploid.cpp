/******************************************************************************/

#include "fbm_byte_acc.h"

/******************************************************************************/

// Estimate nucleotide diversity (pi) at each locus for diploid data. Loci are
// processed in blocks of block_size, checking for user interrupts between
// blocks.

// [[Rcpp::export]]
NumericVector gt_pi_diploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   int ncores,
                                   int block_size = 100000) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  // every individual contributes its dosage and two valid alleles, and all
  // individuals belong to a single group
  std::vector<double> weight(n, 1), copies(n, 2);
  std::vector<size_t> group_ids(n, 0);

  NumericVector pi(m);
  double* p_pi = pi.begin();

  alt_valid_counts(acc, group_ids, 1, weight, copies, ncores, block_size,
                   [&](size_t j, size_t g, double this_allele_count,
                       double this_valid_alleles) {
                     p_pi[j] = (this_valid_alleles > 0) ?
                     (this_allele_count *
                       (this_valid_alleles - this_allele_count) /
                       (this_valid_alleles * (this_valid_alleles -1) /2)) :
                       NA_REAL;
                   });

  return pi;
}

/******************************************************************************/
