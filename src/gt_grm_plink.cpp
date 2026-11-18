/******************************************************************************/

// The Fortran string lengths of BLAS calls (FCONE) need USE_FC_LEN_T before
// the R headers are included
#define USE_FC_LEN_T
#include "pairwise_bitpacked.h"
#include <R_ext/BLAS.h>
// Eigen (from RcppEigen) provides portable, vectorised matrix products that
// do not depend on the BLAS that R uses
#include <Eigen/Core>
#include <chrono>
#include <cmath>
#include <limits>

#ifndef FCONE
#define FCONE
#endif

/******************************************************************************/

// Checked 0-based rows and columns of the FBM for rowInd and colInd
static void grm_offsets(XPtr<FBM> xpBM, const IntegerVector& rowInd,
                        const IntegerVector& colInd,
                        std::vector<std::size_t>& rows,
                        std::vector<std::size_t>& cols) {
  std::size_t nrow_fbm = xpBM->nrow();
  std::size_t ncol_fbm = xpBM->ncol();
  rows.resize(rowInd.size());
  cols.resize(colInd.size());
  for (std::size_t i = 0; i < rows.size(); i++) {
    if (rowInd[i] < 1 || static_cast<std::size_t>(rowInd[i]) > nrow_fbm) {
      Rcpp::stop("Subscript out of bounds.");
    }
    rows[i] = rowInd[i] - 1;
  }
  for (std::size_t l = 0; l < cols.size(); l++) {
    if (colInd[l] < 1 || static_cast<std::size_t>(colInd[l]) > ncol_fbm) {
      Rcpp::stop("Subscript out of bounds.");
    }
    cols[l] = colInd[l] - 1;
  }
}

// Standardised value of each genotype at each locus: 4 values per locus, for
// genotypes 0, 1 and 2 and for missing (always 0). Loci with a zero (or
// undefined) scale, i.e. monomorphic loci, are 0 for every genotype.
// Also fills geno_of_byte: the genotype (0, 1, 2, or 3 for missing) of each
// raw byte of the FBM, decoded through its code256.
static std::vector<double> grm_z_lut(Environment BM,
                                     const NumericVector& center,
                                     const NumericVector& scale,
                                     std::size_t m,
                                     uint8_t* geno_of_byte) {
  if (static_cast<std::size_t>(center.size()) != m ||
      static_cast<std::size_t>(scale.size()) != m) {
    Rcpp::stop("center and scale should have one value per locus.");
  }
  NumericVector code256 = BM["code256"];
  if (code256.size() != 256) {
    Rcpp::stop("The FBM should have a code256 of 256 values.");
  }
  // bytes decoding to anything other than 0, 1 or 2 are missing
  for (int b = 0; b < 256; b++) {
    double g = code256[b];
    geno_of_byte[b] = (g == 0 || g == 1 || g == 2) ?
      static_cast<uint8_t>(g) : 3;
  }
  std::vector<double> zlut(4 * m, 0.0);
  for (std::size_t l = 0; l < m; l++) {
    if (scale[l] > 0 && std::isfinite(scale[l]) && std::isfinite(center[l])) {
      for (int g = 0; g < 3; g++) {
        zlut[4 * l + g] = (g - center[l]) / scale[l];
      }
    }
  }
  return zlut;
}

// Fill the n x nb block Z (column-major) with the standardised genotypes of
// loci l0 to l0 + nb - 1, and return the number of missing genotypes.
// This is done on one thread: the fill is a small part of the work, and an
// OpenMP team spinning between blocks would compete with the threads of the
// BLAS (with 1,000 individuals, half the CPU time went to spinning).
static std::size_t grm_fill_block(double* Z, std::size_t n, std::size_t l0,
                                  std::size_t nb, const unsigned char* p_mat,
                                  std::size_t nrow_fbm,
                                  const std::vector<std::size_t>& rows,
                                  const std::vector<std::size_t>& cols,
                                  const std::vector<double>& zlut,
                                  const uint8_t* geno_of_byte) {
  std::size_t n_missing = 0;
  for (std::size_t s = 0; s < nb; s++) {
    // the bytes of this locus are contiguous in the column-major FBM
    const unsigned char* p_col = p_mat + cols[l0 + s] * nrow_fbm;
    const double* lut = &zlut[4 * (l0 + s)];
    double* z = Z + s * n;
    for (std::size_t i = 0; i < n; i++) {
      uint8_t g = geno_of_byte[p_col[rows[i]]];
      n_missing += (g == 3);
      z[i] = lut[g];
    }
  }
  return n_missing;
}

// G (n x n, lower triangle) += Z Z^T, with Z an n x nb block (column-major).
// With use_blas, the BLAS that R uses (dsyrk): fastest with an optimised
// BLAS (OpenBLAS, MKL, Accelerate), but slow with R's reference BLAS (the
// default on Windows). Otherwise Eigen's own matrix product, in parallel over
// panels of columns of the lower triangle, which needs no BLAS. Eigen uses
// single precision, twice as fast as double, for sub-blocks of at most
// GRM_FLOAT_LOCI loci, whose products are added to the double accumulator:
// the relative error of the relationship matrix is then about 1e-7, well
// below the 6 significant digits that PLINK writes. With the vector
// instructions of every x86-64 CPU (SSE2), this runs at about the speed of
// OpenBLAS in double precision, many times faster than the reference BLAS.
static const std::size_t GRM_FLOAT_LOCI = 1024;

static void grm_rank_update(double* G, const double* Z, std::size_t n,
                            std::size_t nb, bool use_blas, int ncores) {
  if (use_blas) {
    int n_int = static_cast<int>(n);
    int nb_int = static_cast<int>(nb);
    double one = 1.0;
    F77_CALL(dsyrk)("L", "N", &n_int, &nb_int, &one, Z, &n_int, &one, G,
                    &n_int FCONE FCONE);
    return;
  }
  typedef Eigen::Map<Eigen::MatrixXd> MapMat;
  const Eigen::Index n_e = static_cast<Eigen::Index>(n);
  MapMat g(G, n_e, n_e);
  // panels of 256 columns: rows j0 to n - 1 of columns j0 to j0 + w - 1;
  // later panels are shorter, so use a dynamic schedule
  const Eigen::Index panel = 256;
  const Eigen::Index n_panels = (n_e + panel - 1) / panel;
  for (std::size_t s0 = 0; s0 < nb; s0 += GRM_FLOAT_LOCI) {
    const Eigen::Index ns =
      static_cast<Eigen::Index>(std::min(GRM_FLOAT_LOCI, nb - s0));
    // this sub-block of loci in single precision
    Eigen::MatrixXf zf = Eigen::Map<const Eigen::MatrixXd>(
      Z + s0 * n, n_e, ns).cast<float>();
    #pragma omp parallel for num_threads(ncores) schedule(dynamic, 1)
    for (Eigen::Index p = 0; p < n_panels; p++) {
      Eigen::Index j0 = p * panel;
      Eigen::Index w = std::min(panel, n_e - j0);
      Eigen::MatrixXf prod =
        zf.middleRows(j0, n_e - j0) * zf.middleRows(j0, w).transpose();
      g.block(j0, j0, n_e - j0, w) += prod.cast<double>();
    }
  }
}

// Seconds taken by grm_rank_update() for a random n x nb block, to choose
// the faster engine (see gt_pca())
// [[Rcpp::export]]
double grm_rank_update_seconds(int n, int nb, bool use_blas, int ncores) {
  std::vector<double> Z(static_cast<std::size_t>(n) * nb);
  for (std::size_t q = 0; q < Z.size(); q++) {
    Z[q] = static_cast<double>((q * 2654435761u) % 1000) / 500.0 - 1.0;
  }
  std::vector<double> G(static_cast<std::size_t>(n) * n, 0.0);
  auto start = std::chrono::steady_clock::now();
  grm_rank_update(G.data(), Z.data(), n, nb, use_blas, ncores);
  auto end = std::chrono::steady_clock::now();
  return std::chrono::duration<double>(end - start).count();
}

/******************************************************************************/

// Variance-standardised genetic relationship matrix with PLINK 1.9's handling
// of missing genotypes (plink --make-rel, which --pca decomposes):
//   G[i, j] = sum_l z_il z_jl / n_ij,
// where z_il = (g_il - 2 p_l) / sqrt(2 p_l (1 - p_l)) for valid genotypes and
// 0 for missing ones (and at monomorphic loci, which PLINK keeps in n_ij), and
// n_ij is the number of loci valid in both individuals.
// The numerator is accumulated over blocks of block_size loci (lower
// triangle), with the BLAS or with Eigen (use_blas, see grm_rank_update());
// the denominator is counted with the bit-packed popcount kernel of
// pairwise_bitpacked.h.
// center (2 p) and scale (sqrt(2 p (1 - p))) are given for the loci in colInd.
// Returns the full symmetric n x n matrix; pairs with no locus valid in both
// individuals are NaN (as in PLINK). The attribute "n_missing" gives the
// number of missing genotypes.
// [[Rcpp::export]]
NumericMatrix grm_plink_cpp(Environment BM,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd,
                            const NumericVector& center,
                            const NumericVector& scale,
                            int block_size,
                            int ncores,
                            bool use_blas) {
  XPtr<FBM> xpBM = BM["address"];
  const unsigned char* p_mat =
    static_cast<const unsigned char*>(xpBM->matrix());
  std::size_t nrow_fbm = xpBM->nrow();
  std::size_t n = rowInd.size();
  std::size_t m = colInd.size();
  if (block_size < 1) Rcpp::stop("block_size should be positive.");
  std::vector<std::size_t> rows, cols;
  grm_offsets(xpBM, rowInd, colInd, rows, cols);
  uint8_t geno_of_byte[256];
  std::vector<double> zlut = grm_z_lut(BM, center, scale, m, geno_of_byte);

  // numerator: G += Z Z^T for each block Z of standardised genotypes (lower
  // triangle only)
  NumericMatrix G(n, n);
  double* p_G = G.begin();
  std::size_t bs = std::min<std::size_t>(block_size, std::max<std::size_t>(m, 1));
  std::vector<double> Z(n * bs);
  std::size_t n_missing = 0;
  for (std::size_t l0 = 0; l0 < m; l0 += bs) {
    std::size_t nb = std::min(bs, m - l0);
    n_missing += grm_fill_block(Z.data(), n, l0, nb, p_mat, nrow_fbm, rows,
                                cols, zlut, geno_of_byte);
    grm_rank_update(p_G, Z.data(), n, nb, use_blas, ncores);
    Rcpp::checkUserInterrupt();
  }

  // denominator: the number of loci valid in both individuals, for the packed
  // lower triangle (pair (j, i), j >= i, at tri_start[i] + j - i); without
  // missing genotypes, it is the number of loci for every pair
  std::vector<std::size_t> tri_start;
  std::vector<uint32_t> counts;
  if (n_missing > 0) {
    counts = pw_bitpacked_counts<PwValid>(BM, rowInd, colInd, tri_start,
                                          ncores);
  }

  // divide, and copy the lower triangle into the upper one
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = i; j < n; j++) {
      double n_valid = (n_missing > 0) ?
        static_cast<double>(counts[tri_start[i] + j - i]) :
        static_cast<double>(m);
      double value = n_valid > 0 ? p_G[j + i * n] / n_valid : nan;
      p_G[j + i * n] = value;
      p_G[i + j * n] = value;
    }
  }
  // the number of missing genotypes, so that gt_pca() can record whether the
  // genotypes it decomposed were complete (pairwise_grm() drops it)
  G.attr("n_missing") = static_cast<double>(n_missing);
  return G;
}

// Loadings of a gt_pca(): V = Z^T U / d, with Z the standardised
// genotypes (missing as 0) of the loci in colInd, U the n x k eigenvectors
// and d the singular values; for complete data, this is V of the singular
// value decomposition Z = U D V^T.
// [[Rcpp::export]]
NumericMatrix grm_loadings_cpp(Environment BM,
                               const IntegerVector& rowInd,
                               const IntegerVector& colInd,
                               const NumericVector& center,
                               const NumericVector& scale,
                               const NumericMatrix& U,
                               const NumericVector& d,
                               int block_size,
                               int ncores,
                               bool use_blas) {
  XPtr<FBM> xpBM = BM["address"];
  const unsigned char* p_mat =
    static_cast<const unsigned char*>(xpBM->matrix());
  std::size_t nrow_fbm = xpBM->nrow();
  std::size_t n = rowInd.size();
  std::size_t m = colInd.size();
  std::size_t k = U.ncol();
  if (static_cast<std::size_t>(U.nrow()) != n ||
      static_cast<std::size_t>(d.size()) != k) {
    Rcpp::stop("U should be n x k and d of length k.");
  }
  if (block_size < 1) Rcpp::stop("block_size should be positive.");
  std::vector<std::size_t> rows, cols;
  grm_offsets(xpBM, rowInd, colInd, rows, cols);
  uint8_t geno_of_byte[256];
  std::vector<double> zlut = grm_z_lut(BM, center, scale, m, geno_of_byte);

  NumericMatrix V(m, k);
  double* p_V = V.begin();
  const double* p_U = U.begin();
  std::size_t bs = std::min<std::size_t>(block_size, std::max<std::size_t>(m, 1));
  std::vector<double> Z(n * bs);
  std::vector<double> block(bs * k);
  int n_int = static_cast<int>(n);
  int k_int = static_cast<int>(k);
  double one = 1.0, zero = 0.0;
  for (std::size_t l0 = 0; l0 < m; l0 += bs) {
    std::size_t nb = std::min(bs, m - l0);
    grm_fill_block(Z.data(), n, l0, nb, p_mat, nrow_fbm, rows, cols, zlut,
                   geno_of_byte);
    // block (nb x k) = Z^T (nb x n) U (n x k), with the BLAS or Eigen
    int nb_int = static_cast<int>(nb);
    if (use_blas) {
      F77_CALL(dgemm)("T", "N", &nb_int, &k_int, &n_int, &one, Z.data(),
                      &n_int, p_U, &n_int, &zero, block.data(), &nb_int
                      FCONE FCONE);
    } else {
      Eigen::Map<const Eigen::MatrixXd> z(Z.data(), n_int, nb_int);
      Eigen::Map<const Eigen::MatrixXd> u(p_U, n_int, k_int);
      Eigen::Map<Eigen::MatrixXd> out(block.data(), nb_int, k_int);
      out.noalias() = z.transpose() * u;
    }
    // copy into V (m x k, column-major), dividing by d
    for (std::size_t c = 0; c < k; c++) {
      for (std::size_t s = 0; s < nb; s++) {
        p_V[(l0 + s) + c * m] = block[s + c * nb] / d[c];
      }
    }
    Rcpp::checkUserInterrupt();
  }
  return V;
}
