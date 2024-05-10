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
  void increment_king_numerator(Environment k,
                                Environment n_Aa_i,
                                arma::mat& genotype0,
                                arma::mat& genotype1,
                                arma::mat& genotype2,
                                arma::mat& genotype_valid,
                                Environment BM,
                                const IntegerVector& rowInd,
                                const IntegerVector& colInd) {

    arma::mat K = FBM_RW2arma(k);
    arma::mat N_Aa_i = FBM_RW2arma(n_Aa_i);

    XPtr<FBM> xpBM = BM["address"];
    SubBMAcc<unsigned char> macc(xpBM, rowInd, colInd, 1);

  genotype0.zeros();
  genotype1.zeros();
  genotype2.zeros();
  genotype_valid.zeros();

  size_t n = macc.nrow();
  size_t m = macc.ncol();

  for (unsigned int j = 0; j < m; j++) {
    for (unsigned int i = 0; i < n; i++) {
      int value = (macc(i,j));
      if (value == 0){
        genotype0(i, j) = 1;
        genotype_valid(i, j) +=1 ;
      } else if (value==1){
        genotype1(i, j) = 1;
        genotype_valid(i, j) +=1 ;
      } else if (value==2){
        genotype2(i,j) = 1;
        genotype_valid(i, j) +=1 ;
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
  K += ( genotype1 * genotype1.t() -2 *(genotype0 * genotype2.t() + genotype2 * genotype0.t()));
  // use genotype0 to store
  N_Aa_i += (genotype1 * (genotype_valid).t());

}
