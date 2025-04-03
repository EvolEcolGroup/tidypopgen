/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// [[Rcpp::export]]
NumericVector gt_alt_freq_diploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci

  NumericVector freq(m);

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    double this_allele_count = 0; // count of alt alleles at this locus
    double this_valid_alleles = 0; // count of valid alleles at this locus
    for (size_t i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x>-1){
        this_allele_count += x;
        this_valid_alleles +=2;
      }
    }
    freq[j] = (this_valid_alleles > 0) ? (this_allele_count / this_valid_alleles) : NA_REAL;
  }

  return freq;
}

/******************************************************************************/
