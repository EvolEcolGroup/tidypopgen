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

template <class C>
void increment_king_numerator(arma::mat& K,
                          arma::mat& part_temp0,
                          arma::mat& part_temp1,
                          arma::mat& part_temp2,
                          C macc,
                          const IntegerVector& rowInd,
                          const IntegerVector& colInd) {

  // extract the slice of matrix we will operate on
  //part_temp = _extract_submat(macc, part_temp0, part_temp1,part_temp2, rowInd, colInd);

  part_temp0.zeros();
  part_temp1.zeros();
  part_temp2.zeros();

  std::vector<size_t> rows = vec_int_to_size(rowInd, macc.nrow(), 1);
  std::vector<size_t> cols = vec_int_to_size(colInd, macc.ncol(), 1);

  int n = rowInd.size();
  int m = colInd.size();

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      int value = (macc(rows[i], cols[j]));
      if (value == 0){
        part_temp0(i, j) = 1;
      } else if (value==1){
        part_temp1(i, j) = 1;
      } else if (value==2){
        part_temp2(i,j) = 1;
      }
    }}

  // fill the rest with 0s (should be 1 column max)
  int m2 = part_temp0.n_cols;
  if (m2 > m) {
    myassert(m2 == (m + 1), ERROR_BUG);
    for (int i = 0; i < n; i++) {
      part_temp0(i, m) = 0;
      part_temp1(i, m) = 0;
      part_temp2(i, m) = 0;
    }
  }
  K += 2 *( part_temp1 * part_temp1.t() -2 *(part_temp0 * part_temp2.t() + part_temp2 * part_temp0.t()));
  // replace in place to be economical with memory
//  part_temp0 = part_temp0 + part_temp1+part_temp2;
//  K2 += 2* part_temp0 * part_temp0.t();

}

/******************************************************************************/

#define CALL_INCR_KING_NUMERATOR(ACC) {                                                                      \
return increment_king_numerator(armaK, part_temp0, part_temp1, part_temp2, ACC,                      \
                            rowInd, colInd);                                                             \
}                                                                                                        \


// Dispatch function for increment_king_numerator
// [[Rcpp::export]]
void increment_king_numerator(Environment K,
                           arma::mat& part_temp0,
                           arma::mat& part_temp1,
                           arma::mat& part_temp2,
                           Environment BM,
                           const IntegerVector& rowInd,
                           const IntegerVector& colInd) {

  arma::mat armaK = FBM_RW2arma(K);

  DISPATCH_MATACC(CALL_INCR_KING_NUMERATOR)
}

/******************************************************************************/
