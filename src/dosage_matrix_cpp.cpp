#include <bigstatsr/BMCodeAcc.h>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

/******************************************************************************/
// Check a matrix of dosages, and write it into an FBM.code256, used by
// gen_tibble.matrix() and gt_write_fbm_from_dfs().
//
// In R, checking the dosages (with x != round(x), x < 0, x > ploidy and
// is.na(x)) and converting them to raw bytes creates several temporary
// matrices as large as the input (and a double copy of an integer matrix),
// so the peak memory was about 6 times the input matrix. Here the matrix is
// read in place, without any copy.

// Check that the dosages of x (an integer or double matrix) are non-negative
// integers (only if check_integer) no larger than the ploidy of their
// individual (row); ploidy has either one value per row of x or a single
// value for all rows. Missing values (NA or NaN) are allowed. Returns 1 if
// some dosage is not a non-negative integer, otherwise 2 if some dosage is
// larger than its ploidy, and 0 if all dosages are valid (the order of the
// checks of gen_tibble.matrix()).
template <int RTYPE>
int check_dosages(const Matrix<RTYPE>& x, const NumericVector& ploidy,
                  bool check_integer) {
  size_t n = x.nrow();
  size_t m = x.ncol();
  bool one_ploidy = (ploidy.size() == 1);
  if (!one_ploidy && (size_t) ploidy.size() != n) {
    stop("ploidy should have one value per row of x");
  }
  const double* p_ploidy = ploidy.begin();
  typedef typename traits::storage_type<RTYPE>::type stored_t;
  const stored_t* p_x = x.begin();
  bool too_large = false;
  for (size_t j = 0; j < m; j++) {
    const stored_t* x_col = p_x + j * n;
    for (size_t i = 0; i < n; i++) {
      if (traits::is_na<RTYPE>(x_col[i])) continue;
      double value = x_col[i];
      if (check_integer && (value != std::round(value) || value < 0)) {
        return 1;
      }
      if (value > p_ploidy[one_ploidy ? 0 : i]) too_large = true;
    }
  }
  return too_large ? 2 : 0;
}

// [[Rcpp::export]]
int check_dosage_matrix_cpp(SEXP x, const NumericVector& ploidy,
                            bool check_integer) {
  switch (TYPEOF(x)) {
  case INTSXP:
    return check_dosages<INTSXP>(IntegerMatrix(x), ploidy, check_integer);
  case REALSXP:
    return check_dosages<REALSXP>(NumericMatrix(x), ploidy, check_integer);
  default:
    stop("x should be an integer or a numeric matrix");
  }
}

/******************************************************************************/
// Write the dosages of x (an integer or double matrix, with the dimensions of
// the FBM) as bytes into the FBM.code256 BM; missing values are written as
// na_byte. As with as.raw(), values are truncated to integers, and values
// outside 0-255 are written as 0 (dosages are checked beforehand).

template <int RTYPE>
void write_dosages(const Matrix<RTYPE>& x, Environment BM, int na_byte) {
  XPtr<FBM_RW> xpBM = BM["address_rw"];
  if ((size_t) x.nrow() != xpBM->nrow() || (size_t) x.ncol() != xpBM->ncol()) {
    stop("x and the FBM should have the same dimensions");
  }
  if (na_byte < 0 || na_byte > 255) stop("na_byte should be between 0 and 255");
  unsigned char* p_fbm = static_cast<unsigned char*>(xpBM->matrix());
  typedef typename traits::storage_type<RTYPE>::type stored_t;
  const stored_t* p_x = x.begin();
  // both matrices are stored column-major, with the same dimensions
  size_t n_values = (size_t) x.nrow() * (size_t) x.ncol();
  for (size_t k = 0; k < n_values; k++) {
    if (traits::is_na<RTYPE>(p_x[k])) {
      p_fbm[k] = (unsigned char) na_byte;
    } else {
      double value = p_x[k];
      p_fbm[k] = (value >= 0 && value < 256) ? (unsigned char) value : 0;
    }
  }
}

// [[Rcpp::export]]
void write_dosage_matrix_cpp(SEXP x, Environment BM, int na_byte) {
  switch (TYPEOF(x)) {
  case INTSXP:
    write_dosages<INTSXP>(IntegerMatrix(x), BM, na_byte);
    break;
  case REALSXP:
    write_dosages<REALSXP>(NumericMatrix(x), BM, na_byte);
    break;
  default:
    stop("x should be an integer or a numeric matrix");
  }
}

/******************************************************************************/
