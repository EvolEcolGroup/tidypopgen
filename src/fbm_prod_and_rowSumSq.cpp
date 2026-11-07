/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <vector>

/******************************************************************************/

// This is a port of prod_and_rowSumSq from bigsnpr, adapted to work on a
// standard FBM.code256.
//
// For the loci ind (1-based indices into ind_col, center, scale and the rows
// of V), it computes the product of the centred and scaled genotypes (with
// missing genotypes set to 0) by the loadings V, and the sum of squares of
// the scaled genotypes of each individual. The loci are taken in the order
// of ind, so that sums are accumulated in the same order as in the original
// implementation.
//
// If lsq_pcs (1-based column indices of V) is not empty, it also returns,
// for each individual, the sum of v_j v_j' (v_j being row j of
// V[, lsq_pcs]) over the loci at which that individual is missing; this is
// what has to be removed from crossprod(V[, lsq_pcs]) to obtain the
// crossproduct over the loci genotyped in that individual, as needed by the
// least squares projection.
//
// For speed, genotypes are decoded with a local copy of code256, the
// loadings of each locus are copied into a small contiguous vector, and the
// products are accumulated in an individual-major (K x n) array, so that the
// innermost loop runs over contiguous memory.

// [[Rcpp::export]]
List fbm256_prod_and_rowSumsSq(Environment BM,
                               const IntegerVector& ind_row,
                               const IntegerVector& ind_col,
                               const NumericVector& center,
                               const NumericVector& scale,
                               const NumericMatrix& V,
                               const IntegerVector& ind,
                               const IntegerVector& lsq_pcs) {

  XPtr<FBM> xpBM = BM["address"];
  // raw (undecoded) bytes of the FBM, stored column-major
  const unsigned char* p_mat =
    static_cast<const unsigned char*>(xpBM->matrix());
  std::size_t nrow_fbm = xpBM->nrow();
  std::size_t ncol_fbm = xpBM->ncol();

  // local copy of the code to decode the raw bytes
  NumericVector code256 = BM["code256"];
  double code[256];
  for (int b = 0; b < 256; b++) code[b] = code256[b];

  std::size_t n = ind_row.size(); // number of individuals
  std::size_t m = ind.size(); // number of loci used in this call
  std::size_t n_v = V.nrow(); // number of loci in V (and ind_col)
  std::size_t K = V.ncol(); // number of components
  std::size_t K2 = lsq_pcs.size(); // number of components for least squares
  myassert_size(ind_col.size(), n_v);
  myassert_size(center.size(), n_v);
  myassert_size(scale.size(), n_v);

  // 0-based offsets of the rows (checked against the dimensions of the FBM)
  std::vector<std::size_t> rows(n);
  for (std::size_t i = 0; i < n; i++) {
    if (ind_row[i] < 1 || static_cast<std::size_t>(ind_row[i]) > nrow_fbm) {
      Rcpp::stop("Subscript out of bounds.");
    }
    rows[i] = ind_row[i] - 1;
  }
  // 0-based columns of V used for least squares
  std::vector<std::size_t> pcs(K2);
  for (std::size_t k = 0; k < K2; k++) {
    if (lsq_pcs[k] < 1 || static_cast<std::size_t>(lsq_pcs[k]) > K) {
      Rcpp::stop("lsq_pcs out of bounds.");
    }
    pcs[k] = lsq_pcs[k] - 1;
  }

  const double* p_V = V.begin();
  // accumulators, individual-major: K (or K2 * K2) values per individual
  std::vector<double> xv_t(n * K, 0.0);
  std::vector<double> row_sums_sq(n, 0.0);
  std::vector<double> miss_t(n * K2 * K2, 0.0);
  // loadings of the current locus (all components, and those for lsq)
  std::vector<double> v_j(K), v_j_lsq(K2);

  for (std::size_t jj = 0; jj < m; jj++) {
    // index of this locus in ind_col, center, scale and V
    if (ind[jj] < 1 || static_cast<std::size_t>(ind[jj]) > n_v) {
      Rcpp::stop("Subscript out of bounds.");
    }
    std::size_t j = ind[jj] - 1;
    if (ind_col[j] < 1 || static_cast<std::size_t>(ind_col[j]) > ncol_fbm) {
      Rcpp::stop("Subscript out of bounds.");
    }
    const unsigned char* p_col = p_mat + (ind_col[j] - 1) * nrow_fbm;
    double center_j = center[j];
    double scale_j = scale[j];
    for (std::size_t k = 0; k < K; k++) v_j[k] = p_V[j + k * n_v];
    for (std::size_t k = 0; k < K2; k++) v_j_lsq[k] = v_j[pcs[k]];

    for (std::size_t i = 0; i < n; i++) {
      double x = code[p_col[rows[i]]];
      if (x > -1) {
        // here we need to center and standardise
        x = (x - center_j) / scale_j;
      } else {
        // missing genotypes are set to 0 (i.e. the mean)
        x = 0;
        // and, for least squares, their loadings are recorded (miss_t is
        // empty without least squares, so it must not be indexed then)
        if (K2 > 0) {
          double* p_miss = miss_t.data() + i * K2 * K2;
          for (std::size_t k1 = 0; k1 < K2; k1++) {
            for (std::size_t k2 = 0; k2 < K2; k2++) {
              p_miss[k1 + k2 * K2] += v_j_lsq[k1] * v_j_lsq[k2];
            }
          }
        }
      }
      row_sums_sq[i] += x * x;
      double* p_xv = xv_t.data() + i * K;
      for (std::size_t k = 0; k < K; k++) {
        p_xv[k] += x * v_j[k];
      }
    }
  }

  // copy the accumulators into individual-by-column matrices
  NumericMatrix XV(n, K);
  NumericMatrix miss(n, K2 * K2);
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t k = 0; k < K; k++) XV(i, k) = xv_t[i * K + k];
    for (std::size_t k = 0; k < K2 * K2; k++) {
      miss(i, k) = miss_t[i * K2 * K2 + k];
    }
  }

  return List::create(XV, NumericVector(row_sums_sq.begin(),
                                        row_sums_sq.end()), miss);
}

/******************************************************************************/
