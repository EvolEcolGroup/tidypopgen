/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// Grouped (per-population) version of gt_pi_polyploid(); see that function
// for details of the estimator. Kept as a separate function from
// gt_grouped_pi_diploid() so the diploid fast path stays free of the extra
// ploidy lookup this requires.

// [[Rcpp::export]]
ListOf<NumericMatrix> gt_grouped_pi_polyploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   int ngroups,
                                   const NumericVector& ploidy,
                                   int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci

  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }

  NumericMatrix pi(m, ngroups);
  NumericMatrix valid_alleles(m, ngroups);

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x > -1){
        pi(j, groupIds[i]) += x;
        valid_alleles(j, groupIds[i]) += ploidy[i];
      }
    }
    // now for each group, divide freq by valid_alleles
    for (int group_i = 0; group_i < ngroups; group_i++) {
      double this_valid_alleles = valid_alleles(j, group_i);
      // guard against 0 or 1 valid allele copies in this group at this
      // locus, where the denominator would be 0
      pi(j, group_i) = (this_valid_alleles > 1) ?
        (pi(j, group_i) * (this_valid_alleles - pi(j, group_i)) /
          (this_valid_alleles * (this_valid_alleles - 1) / 2)) : NA_REAL;
    }
  }

  return List::create(_["pi"]  = pi,
                      _["n"] = valid_alleles);
}

/******************************************************************************/
