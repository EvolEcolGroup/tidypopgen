/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// [[Rcpp::export]]
NumericVector cpp_alt_freq_dip_pseudo(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const NumericVector& ploidy,
                                   const bool is_pseudohap,
                                   int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci

  // multiplier to convert dosages for pseudohaploid data
  NumericVector pseudo_multiplier (n);
  for (size_t i = 0; i < n; i++) {
    if (is_pseudohap) {
      pseudo_multiplier[i] = 1 / (3 - ploidy[i]);
    } else {
      pseudo_multiplier[i] = 1;
    }
  }

  NumericVector freq(m);

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) { // loop over loci
    double this_allele_count = 0; // count of alt alleles at this locus
    double this_valid_alleles = 0; // count of valid alleles at this locus
    for (size_t i = 0; i < n; i++) { // loop over individuals
      double x = macc(i, j);
      if (x > -1){
        this_allele_count += x * pseudo_multiplier[i];
        this_valid_alleles += ploidy(i);
      }
    }
    freq[j] = (this_valid_alleles > 0) ? (this_allele_count / this_valid_alleles) : NA_REAL;
  }

  return freq;
}

/******************************************************************************/
