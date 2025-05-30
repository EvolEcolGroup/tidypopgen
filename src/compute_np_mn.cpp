#include <Rcpp.h>
using namespace Rcpp;

// function to compute the number of non-missing values (np) and the mean of 
// the inverse (mn) for each row in a matrix. This is used as part of loci_globals_stats

// [[Rcpp::export]]
List compute_np_mn(NumericMatrix n) {
  int nrow = n.nrow();
  int ncol = n.ncol();
  
  NumericVector np(nrow, 0.0);
  NumericVector denom(nrow, 0.0);
  
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      double val = n(i, j);
      if (!NumericVector::is_na(val)) {
        np[i] += 1.0;
        denom[i] += 1.0 / val;
      }
    }
  }
  
  NumericVector mn(nrow);
  for (int i = 0; i < nrow; ++i) {
    mn[i] = (denom[i] > 0.0) ? np[i] / denom[i] : NA_REAL;
  }
  
  return List::create(Named("np") = np,
                      Named("mn") = mn);
}
