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
  void increment_ibs_counts(Environment k,
                            Environment k2,
                            arma::mat& genotype0,
                            arma::mat& genotype1,
                            arma::mat& genotype2,
                            Environment BM,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd) {

    arma::mat K = FBM_RW2arma(k);
    arma::mat K2 = FBM_RW2arma(k2);

  XPtr<FBM> xpBM = BM["address"];
  SubBMAcc<unsigned char> macc(xpBM, rowInd, colInd, 1);


  genotype0.zeros();
  genotype1.zeros();
  genotype2.zeros();

  size_t n = macc.nrow();
  size_t m = macc.ncol();

  for (unsigned int j = 0; j < m; j++) {
    for (unsigned int i = 0; i < n; i++) {
      int value = (macc(i,j));
      if (value == 0){
        genotype0(i, j) = 1;
      } else if (value==1){
        genotype1(i, j) = 1;
      } else if (value==2){
        genotype2(i,j) = 1;
      }
    }}

  // fill the rest with 0s (should be 1 column max)
  size_t m2 = genotype0.n_cols;
  if (m2 > m) {
    myassert(m2 == (m + 1), ERROR_BUG);
    for (unsigned int i = 0; i < n; i++) {
      genotype0(i, m) = 0;
      genotype1(i, m) = 0;
      genotype2(i, m) = 0;
    }
  }
  K += 2 *( genotype2 * genotype2.t() + genotype1 * genotype1.t() + genotype0 * genotype0.t()) +
  genotype1 * (genotype0 + genotype2).t() + (genotype0 + genotype2) * genotype1.t();
  // replace in place to be economical with memory
  // this is to count the available genotypes (i.e. non-missing)
  genotype0 = genotype0 + genotype1 + genotype2;
  K2 += 2* genotype0 * genotype0.t();

}

/******************************************************************************/

