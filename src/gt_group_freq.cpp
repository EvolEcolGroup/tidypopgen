/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// [[Rcpp::export]]
ListOf<NumericMatrix> gt_grouped_alt_freq_diploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   int ngroups,
                                   int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci

  NumericMatrix freq(ngroups, m);
  NumericMatrix valid_alleles(ngroups, m);

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x>-1){
        freq(groupIds[i],j) += x;
        valid_alleles(groupIds[i],j) +=2;
      }
    }
    // now for each group, divide freq by valid_alleles
    for (size_t group_i = 0; group_i < ngroups; group_i++) {
      freq(group_i,j) = freq(group_i,j) / valid_alleles(group_i,j);
    }
  }

  return List::create(_["freq_alt"]  = freq,
                      _["n"] = valid_alleles);
}

/******************************************************************************/
