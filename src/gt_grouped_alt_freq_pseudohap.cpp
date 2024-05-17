/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// [[Rcpp::export]]
ListOf<NumericMatrix> gt_grouped_alt_freq_pseudohap(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   int ngroups,
                                   const IntegerVector& ploidy,
                                   int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci

  NumericMatrix freq(m, ngroups);
  NumericMatrix valid_alleles(m, ngroups);

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x>-1){
        if (ploidy[i]==2){
          freq(j, groupIds[i]) += x;
          valid_alleles(j, groupIds[i]) +=2;
        } else {
          freq(j, groupIds[i]) += x/2;
          valid_alleles(j, groupIds[i]) +=1;
        }
      }
    }
    // now for each group, divide freq by valid_alleles
    for (size_t group_i = 0; group_i < ngroups; group_i++) {
      freq(j, group_i) = freq(j, group_i) / valid_alleles(j, group_i);
    }
  }

  return List::create(_["freq_alt"]  = freq,
                      _["n"] = valid_alleles);
}

/******************************************************************************/
