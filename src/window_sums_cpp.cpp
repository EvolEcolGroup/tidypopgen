#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

/******************************************************************************/
// Sums of the non-missing values in each window, for each column of x.
//
// This replaces runner::sum_run() and runner::mean_run() (with na_rm = TRUE)
// in windows_stats_generic(): the windows are defined once (see
// .window_definitions()), and each window is then summed directly, which is
// much faster than finding the loci of each window with runner (especially
// when windows are defined in bp), and works on many columns at once.
//
// x: a matrix with loci as rows (in the order of the loci table)
// loci: the 1-based indices of the loci (rows of x), ordered by chromosome
// lo, hi: window w includes loci[lo[w]], ..., loci[hi[w] - 1] (0-based)
// long_double: should sums be accumulated in long double (as R's sum()),
//   rather than in double (as runner::sum_run())?
// Returns a list with two matrices with one row per window and one column per
// column of x: `sum`, the sum of the non-missing values (accumulated in the
// order of the loci), and `n`, the number of non-missing values.

// [[Rcpp::export]]
List window_sums_cpp(const NumericMatrix& x,
                     const IntegerVector& loci,
                     const IntegerVector& lo,
                     const IntegerVector& hi,
                     bool long_double) {
  size_t n_loci = x.nrow();
  size_t n_col = x.ncol();
  size_t n_win = lo.size();
  if ((size_t) hi.size() != n_win) stop("lo and hi should have the same length");
  // check the indices once, so that we can then use raw pointers
  for (R_xlen_t k = 0; k < loci.size(); k++) {
    if (loci[k] < 1 || (size_t) loci[k] > n_loci) {
      stop("loci should be indices of rows of x");
    }
  }
  for (size_t w = 0; w < n_win; w++) {
    if (lo[w] < 0 || hi[w] < lo[w] || hi[w] > loci.size()) {
      stop("invalid window boundaries");
    }
  }
  NumericMatrix sums(n_win, n_col);
  NumericMatrix counts(n_win, n_col);
  const double* p_x = x.begin();
  const int* p_loci = loci.begin();

  for (size_t s = 0; s < n_col; s++) {
    const double* x_col = p_x + s * n_loci;
    for (size_t w = 0; w < n_win; w++) {
      double count = 0;
      if (long_double) {
        long double sum = 0.0;
        for (int k = lo[w]; k < hi[w]; k++) {
          double value = x_col[p_loci[k] - 1];
          if (!std::isnan(value)) {
            sum += value;
            count++;
          }
        }
        sums(w, s) = (double) sum;
      } else {
        double sum = 0.0;
        for (int k = lo[w]; k < hi[w]; k++) {
          double value = x_col[p_loci[k] - 1];
          if (!std::isnan(value)) {
            sum += value;
            count++;
          }
        }
        sums(w, s) = sum;
      }
      counts(w, s) = count;
    }
  }
  return List::create(_["sum"] = sums, _["n"] = counts);
}

/******************************************************************************/
