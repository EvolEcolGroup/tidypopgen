/******************************************************************************/

#include "fbm_byte_acc.h"

/******************************************************************************/

// Estimate nucleotide diversity (pi) at each locus, generalised to
// arbitrary and mixed ploidy. Mathematically this is the unbiased gene
// diversity estimator of Nei & Roychoudhury (1974) / DeGiorgio, Jankovic &
// Rosenberg (2010) eqn 3: M/(M-1) * (1 - sum_i p_i^2), where M is the total
// number of allele copies. When the dataset is marked pseudohaploid,
// individuals with ploidy one retain 0/2 storage while diploid individuals
// retain their usual dosage. The dataset-level flag distinguishes this from
// ordinary 0/1 haploid calls. Individuals of different ploidy simply use
// different lookup tables (see fbm_byte_acc.h), so this is as fast as the
// diploid version.

// [[Rcpp::export]]
NumericVector gt_pi_polyploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const NumericVector& ploidy,
                                   bool is_pseudohaploid,
                                   int ncores,
                                   int block_size = 100000) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }

  // pseudohaploid individuals contribute half their (0/2) dosage, and every
  // individual contributes as many valid alleles as its ploidy
  std::vector<double> weight(n), copies(n);
  for (size_t i = 0; i < n; i++) {
    weight[i] = (is_pseudohaploid && ploidy[i] == 1) ? 0.5 : 1;
    copies[i] = ploidy[i];
  }
  // all individuals belong to a single group
  std::vector<size_t> group_ids(n, 0);

  NumericVector pi(m);
  double* p_pi = pi.begin();

  alt_valid_counts(acc, group_ids, 1, weight, copies, ncores, block_size,
                   [&](size_t j, size_t g, double this_allele_count,
                       double this_valid_alleles) {
                     // guard against a single valid allele copy, where
                     // valid-1 would be 0
                     p_pi[j] = (this_valid_alleles > 1) ?
                     (this_allele_count *
                       (this_valid_alleles - this_allele_count) /
                       (this_valid_alleles * (this_valid_alleles -1) /2)) :
                       NA_REAL;
                   });

  return pi;
}

/******************************************************************************/
