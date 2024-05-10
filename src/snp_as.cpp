/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/BMAcc-dispatcher.h>

/******************************************************************************/

inline arma::mat FBM_RW2arma(Rcpp::Environment BM) {

  Rcpp::XPtr<FBM_RW> xpBM = BM["address_rw"];
  myassert(xpBM->matrix_type() == 8,
           "Mapping to arma::mat is available for 'double' FBMs only.");

  return arma::mat((double*)xpBM->matrix(), xpBM->nrow(), xpBM->ncol(), false);
}

/******************************************************************************/

/******************************************************************************/

  // [[Rcpp::export]]
  void increment_as_counts(Environment k, // to store tcrossprod(dos_mat-1)
                            Environment k2, // to store tcrossprod(na_mat)
                            arma::mat& na_mat,
                            arma::mat& dos_mat,
                            Environment BM,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd) {

    arma::mat K = FBM_RW2arma(k);
    arma::mat K2 = FBM_RW2arma(k2);

  XPtr<FBM> xpBM = BM["address"];
  SubBMAcc<unsigned char> macc(xpBM, rowInd, colInd, 1);

  dos_mat.zeros();
  na_mat.zeros();
  na_mat = na_mat +1; // just reset it to one dirctly!!!


  size_t n = macc.nrow();
  size_t m = macc.ncol();

  for (unsigned int j = 0; j < m; j++) {
    for (unsigned int i = 0; i < n; i++) {
      int value = (macc(i,j));
      if (value <3){
        dos_mat(i, j) = value-1; // note that we want the matrix of (dos-1)
      } else {
        dos_mat(i, j) = 0;
        na_mat(i, j) = 0;
      }
    }}

  // fill the rest with 0s (should be 1 column max)
  size_t m2 = dos_mat.n_cols;
  if (m2 > m) {
    myassert(m2 == (m + 1), ERROR_BUG);
    for (unsigned int i = 0; i < n; i++) {
      dos_mat(i, m) = 1;
      na_mat(i, m) = 0;
    }
  }
  K += dos_mat * dos_mat.t();
  K2 += na_mat * na_mat.t();

}

/******************************************************************************/

