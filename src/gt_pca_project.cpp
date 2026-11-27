/******************************************************************************/

// The Fortran string lengths of BLAS calls (FCONE) need USE_FC_LEN_T before
// the R headers are included
#define USE_FC_LEN_T
#include <bigstatsr/BMAcc.h>
#include <R_ext/BLAS.h>
// Eigen (from RcppEigen) provides matrix products that do not depend on the
// BLAS that R uses (see gt_grm_plink.cpp)
#include <Eigen/Core>
#include <cmath>
#include <vector>

#ifndef FCONE
#define FCONE
#endif

/******************************************************************************/

// Projection of individuals on the loadings V of a gt_pca, for
// predict.gt_pca().
//
// For each individual, it computes over the loci of colInd:
// - XV: the product of the centred and scaled genotypes (missing genotypes
//   set to 0) by the loadings V (n x K);
// - the sum of squares of the scaled genotypes (for OADP);
// - the number of missing genotypes (for the "plink" projection, which
//   divides XV by the fraction of loci genotyped);
// - if lsq_pcs (1-based columns of V) is not empty, the sum of v_l v_l' over
//   the loci l missing in that individual, with v_l row l of V[, lsq_pcs]
//   (n x K2^2, as vectorised K2 x K2 matrices): what has to be removed from
//   crossprod(V[, lsq_pcs]) to obtain the crossproduct over the loci
//   genotyped in that individual, for the least squares projection.
//
// Rather than visiting each genotype in turn, the loci are processed in
// blocks: the standardised genotypes Z (n x nb) and the indicators of missing
// genotypes M (n x nb) of a block are decoded, and XV += Z V_block and
// Miss += M W_block are matrix products, with W_block (nb x K2(K2+1)/2) the
// products of the loadings of each locus. The products use the BLAS
// (use_blas) or Eigen, as for gt_pca(). V has one row per locus of colInd,
// with center and scale (from the gt_pca object).
// [[Rcpp::export]]
List pca_project_cpp(Environment BM,
                     const IntegerVector& rowInd,
                     const IntegerVector& colInd,
                     const NumericVector& center,
                     const NumericVector& scale,
                     const NumericMatrix& V,
                     const IntegerVector& lsq_pcs,
                     int block_size,
                     bool use_blas,
                     int ncores) {
  XPtr<FBM> xpBM = BM["address"];
  // raw (undecoded) bytes of the FBM, stored column-major
  const unsigned char* p_mat =
    static_cast<const unsigned char*>(xpBM->matrix());
  const std::size_t nrow_fbm = xpBM->nrow();
  const std::size_t ncol_fbm = xpBM->ncol();
  const std::size_t n = rowInd.size();
  const std::size_t m = colInd.size();
  const std::size_t K = V.ncol();
  const std::size_t K2 = lsq_pcs.size();
  // the distinct products of the loadings of a locus (k1 <= k2)
  const std::size_t P = K2 * (K2 + 1) / 2;
  if (static_cast<std::size_t>(V.nrow()) != m ||
      static_cast<std::size_t>(center.size()) != m ||
      static_cast<std::size_t>(scale.size()) != m) {
    Rcpp::stop("V, center and scale should have one row or value per locus.");
  }
  if (block_size < 1) Rcpp::stop("block_size should be positive.");

  // 0-based rows, columns and least squares components, checked
  std::vector<std::size_t> rows(n), cols(m), pcs(K2);
  for (std::size_t i = 0; i < n; i++) {
    if (rowInd[i] < 1 || static_cast<std::size_t>(rowInd[i]) > nrow_fbm) {
      Rcpp::stop("Subscript out of bounds.");
    }
    rows[i] = rowInd[i] - 1;
  }
  for (std::size_t l = 0; l < m; l++) {
    if (colInd[l] < 1 || static_cast<std::size_t>(colInd[l]) > ncol_fbm) {
      Rcpp::stop("Subscript out of bounds.");
    }
    cols[l] = colInd[l] - 1;
  }
  for (std::size_t k = 0; k < K2; k++) {
    if (lsq_pcs[k] < 1 || static_cast<std::size_t>(lsq_pcs[k]) > K) {
      Rcpp::stop("lsq_pcs out of bounds.");
    }
    pcs[k] = lsq_pcs[k] - 1;
  }

  // decoded value of each raw byte (NA, or anything below -1, is missing)
  NumericVector code256 = BM["code256"];
  if (code256.size() != 256) {
    Rcpp::stop("The FBM should have a code256 of 256 values.");
  }
  double code[256];
  bool missing_byte[256];
  for (int b = 0; b < 256; b++) {
    code[b] = code256[b];
    missing_byte[b] = !(code[b] > -1);
  }

  const double* p_V = V.begin();
  const std::size_t bs =
    std::min<std::size_t>(block_size, std::max<std::size_t>(m, 1));
  std::vector<double> Z(n * bs), Vb(bs * K);
  std::vector<double> Mi(K2 > 0 ? n * bs : 0), Wb(K2 > 0 ? bs * P : 0);
  NumericMatrix XV(n, K);
  double* p_XV = XV.begin();
  std::vector<double> miss_u(n * P, 0.0);
  std::vector<double> row_sums_sq(n, 0.0);
  std::vector<double> n_missing(n, 0.0);

  for (std::size_t l0 = 0; l0 < m; l0 += bs) {
    const std::size_t nb = std::min(bs, m - l0);
    // decode the block (one thread: see grm_fill_block() in gt_grm_plink.cpp)
    for (std::size_t s = 0; s < nb; s++) {
      const std::size_t l = l0 + s;
      const unsigned char* p_col = p_mat + cols[l] * nrow_fbm;
      // loci without a scale (monomorphic) are 0 for every genotype
      const bool scaled = scale[l] > 0 && std::isfinite(scale[l]);
      const double c = center[l];
      const double inv_s = scaled ? 1.0 / scale[l] : 0.0;
      double* z = &Z[s * n];
      double* mi = (K2 > 0) ? &Mi[s * n] : nullptr;
      for (std::size_t i = 0; i < n; i++) {
        const unsigned char byte = p_col[rows[i]];
        const bool is_missing = missing_byte[byte];
        const double v = is_missing ? 0.0 : (code[byte] - c) * inv_s;
        z[i] = v;
        row_sums_sq[i] += v * v;
        n_missing[i] += is_missing;
        if (mi) mi[i] = is_missing;
      }
      // the loadings of this locus, and the products of those for lsq
      for (std::size_t k = 0; k < K; k++) Vb[s + k * nb] = p_V[l + k * m];
      if (K2 > 0) {
        std::size_t q = 0;
        for (std::size_t a = 0; a < K2; a++) {
          for (std::size_t b = a; b < K2; b++, q++) {
            Wb[s + q * nb] = p_V[l + pcs[a] * m] * p_V[l + pcs[b] * m];
          }
        }
      }
    }
    // XV += Z V_block, and Miss += M W_block
    if (use_blas) {
      int n_i = static_cast<int>(n), nb_i = static_cast<int>(nb);
      int K_i = static_cast<int>(K), P_i = static_cast<int>(P);
      double one = 1.0;
      F77_CALL(dgemm)("N", "N", &n_i, &K_i, &nb_i, &one, Z.data(), &n_i,
                      Vb.data(), &nb_i, &one, p_XV, &n_i FCONE FCONE);
      if (K2 > 0) {
        F77_CALL(dgemm)("N", "N", &n_i, &P_i, &nb_i, &one, Mi.data(), &n_i,
                        Wb.data(), &nb_i, &one, miss_u.data(), &n_i
                        FCONE FCONE);
      }
    } else {
      typedef Eigen::Map<Eigen::MatrixXd> MapMat;
      typedef Eigen::Map<const Eigen::MatrixXd> ConstMapMat;
      const Eigen::Index n_e = static_cast<Eigen::Index>(n);
      const Eigen::Index nb_e = static_cast<Eigen::Index>(nb);
      MapMat xv(p_XV, n_e, static_cast<Eigen::Index>(K));
      MapMat mu(miss_u.data(), n_e, static_cast<Eigen::Index>(P));
      ConstMapMat z(Z.data(), n_e, nb_e);
      ConstMapMat vb(Vb.data(), nb_e, static_cast<Eigen::Index>(K));
      ConstMapMat mi(Mi.data(), K2 > 0 ? n_e : 0, K2 > 0 ? nb_e : 0);
      ConstMapMat wb(Wb.data(), K2 > 0 ? nb_e : 0,
                     static_cast<Eigen::Index>(P));
      // in parallel over blocks of individuals (rows)
      const Eigen::Index chunk = 256;
      const Eigen::Index n_chunks = (n_e + chunk - 1) / chunk;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic, 1)
      for (Eigen::Index c = 0; c < n_chunks; c++) {
        Eigen::Index i0 = c * chunk;
        Eigen::Index ni = std::min(chunk, n_e - i0);
        xv.middleRows(i0, ni).noalias() += z.middleRows(i0, ni) * vb;
        if (K2 > 0) {
          mu.middleRows(i0, ni).noalias() += mi.middleRows(i0, ni) * wb;
        }
      }
    }
    Rcpp::checkUserInterrupt();
  }

  // the missing crossproducts as full (vectorised) K2 x K2 matrices
  NumericMatrix miss(n, K2 * K2);
  std::size_t q = 0;
  for (std::size_t a = 0; a < K2; a++) {
    for (std::size_t b = a; b < K2; b++, q++) {
      for (std::size_t i = 0; i < n; i++) {
        miss(i, a + b * K2) = miss_u[i + q * n];
        miss(i, b + a * K2) = miss_u[i + q * n];
      }
    }
  }
  return List::create(
    _["XV"] = XV,
    _["row_sums_sq"] = NumericVector(row_sums_sq.begin(), row_sums_sq.end()),
    _["n_missing"] = NumericVector(n_missing.begin(), n_missing.end()),
    _["miss"] = miss
  );
}
