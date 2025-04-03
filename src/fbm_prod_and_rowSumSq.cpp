/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// This is a port of prod_and_rowSumSq from bigsnpr, adapted to work on a standard FBM

// [[Rcpp::export]]
List fbm256_prod_and_rowSumsSq(Environment BM,
                        const IntegerVector& ind_row,
                        const IntegerVector& ind_col,
                        const NumericVector& center,
                        const NumericVector& scale,
                        const NumericMatrix& V) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, ind_row, ind_col, BM["code256"], 1);


  size_t n = macc.nrow(); //Number of individuals
  size_t m = macc.ncol(); //Number of sites
  myassert_size(m, V.rows()); //check number of sites same as number in V from PCA
  size_t K = V.cols();
  size_t i, j, k;

  NumericMatrix XV(n, K);
  NumericVector rowSumsSq(n);

  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      double x = macc(i, j); // here we need to center and standardise
      if (x > -1){
        x = (x-center[j])/scale[j];
      } else {
//        Rcout<<"impute"<<std::endl;
        x = 0;
      }
      rowSumsSq[i] += x*x;
      for (k = 0; k < K; k++) {
        XV(i, k) += x * V(j, k);
      }
    }
  }

  return List::create(XV, rowSumsSq);
}

