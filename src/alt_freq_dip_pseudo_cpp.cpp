/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix alt_freq_dip_pseudo_cpp(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const NumericVector& ploidy,
                                   int ncores,
                                   bool as_counts) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci

  // multiplier to convert dosages for pseudohaploid data
  NumericVector pseudo_multiplier (n);
  for (size_t i = 0; i < n; i++) {
    pseudo_multiplier[i] = 1 / (3 - ploidy[i]);
  }

  // Matrix to store frequency and valid alleles in two columns
  // the first column initially stores counts of alternate alleles and
  // can be returned if return_counts = TRUE
  NumericMatrix freq_mat(m,2);

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) { // loop over loci
    for (size_t i = 0; i < n; i++) { // loop over individuals
      double x = macc(i, j);
      if (x > -1){
        freq_mat(j,0) += x * pseudo_multiplier[i];
        freq_mat(j,1)  += ploidy(i);
      }
    }
  }
  
  if (as_counts){
    colnames(freq_mat) = CharacterVector::create("n_alt", "n_valid");
    return freq_mat;
  }
  
  for (size_t j = 0; j < m; j++) { // loop over loci
    if (freq_mat(j,1)>0){
      freq_mat(j,0) = freq_mat(j,0) / freq_mat(j,1);
    } else {
      freq_mat(j,0) = NA_REAL;
    }
  }
  
  colnames(freq_mat) = CharacterVector::create("freq", "n_valid");
  return freq_mat;
}

/******************************************************************************/
