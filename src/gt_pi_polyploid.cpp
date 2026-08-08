/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// Estimate nucleotide diversity (pi) at each locus, generalised to
// arbitrary and mixed ploidy. Mathematically this is the unbiased gene
// diversity estimator of Nei & Roychoudhury (1974) / DeGiorgio, Jankovic &
// Rosenberg (2010) eqn 3: M/(M-1) * (1 - sum_i p_i^2), where M is the total
// number of allele copies. Kept as a separate function from gt_pi_diploid()
// so the diploid fast path stays free of the extra ploidy lookup this
// requires.

// [[Rcpp::export]]
NumericVector gt_pi_polyploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const NumericVector& ploidy,
                                   int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci

  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }

  NumericVector pi(m);

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    double this_allele_count = 0; // count of alt alleles at this locus
    double this_valid_alleles = 0; // count of valid alleles at this locus
    for (size_t i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x > -1){
        this_allele_count += x;
        this_valid_alleles += ploidy[i];
      }
    }
    // guard against a single valid allele copy, where valid-1 would be 0
    pi[j] = (this_valid_alleles > 1) ?
    (this_allele_count * (this_valid_alleles - this_allele_count) /
      (this_valid_alleles * (this_valid_alleles -1) /2)) : NA_REAL;
  }

  return pi;
}

/******************************************************************************/
