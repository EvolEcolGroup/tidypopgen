/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// [[Rcpp::export]]
ListOf<NumericMatrix> gt_grouped_pi_diploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   int ngroups,
                                   int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci

  NumericMatrix pi(m, ngroups);
  NumericMatrix valid_alleles(m, ngroups);

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x > -1){
        pi(j, groupIds[i]) += x;
        valid_alleles(j, groupIds[i]) +=2;
      }
    }
    // now for each group, divide freq by valid_alleles
    for (int group_i = 0; group_i < ngroups; group_i++) {
      pi(j, group_i) = (pi(j, group_i) * (valid_alleles(j, group_i) - pi(j, group_i)) / 
        (valid_alleles(j, group_i) * (valid_alleles(j, group_i) -1) /2));
    }
  }

  return List::create(_["pi"]  = pi,
                      _["n"] = valid_alleles);
}

/******************************************************************************/
