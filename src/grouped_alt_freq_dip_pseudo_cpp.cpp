/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix grouped_alt_freq_dip_pseudo_cpp(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   int ngroups,
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
  
  // We use a single large matrix, the first ngroup columns are frequencies
  // and the next ngroups columns are valid alleles. If as_counts = TRUE, then
  // we return the counts before we compute the frequencies
  
  NumericMatrix freq_mat(m, ngroups*2);

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) { // loop over loci
    for (size_t i = 0; i < n; i++) { //loop over individuals
      double x = macc(i, j);
      if (x > -1){
        freq_mat(j, groupIds[i]) += x * pseudo_multiplier[i];
        freq_mat(j, ngroups+groupIds[i]) += ploidy[i];
      }
    }
  }
  
  if (as_counts){
    return freq_mat;
  }
  // compute frequencies
  for (size_t j = 0; j < m; j++) { // loop over loci
    for (int group_i = 0; group_i < ngroups; group_i++) {
      freq_mat(j, group_i) = freq_mat(j, group_i) / 
        freq_mat(j, ngroups + group_i); // divide alt count by valid count
    }
  }
  
  return freq_mat;
}

/******************************************************************************/
