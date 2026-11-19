/******************************************************************************/

// Linkage disequilibrium (r^2) between pairs of loci with missing genotypes,
// for loci_ld_clump() (bigsnpr's clumping).
//
// r^2 is the squared correlation of the allele counts (0, 1, 2) of two loci
// over the individuals genotyped at both, as in PLINK 1.9 (--r2,
// --indep-pairwise); with complete genotypes, it is the r^2 of bigsnpr's
// snp_clumping(). Each locus is packed into three bit planes over the
// individuals (64 per word): genotyped, genotype >= 1, and genotype == 2, so
// that the genotype is the sum of the last two planes and its square is the
// first of them plus 3 times the second. For a pair of loci a and b, the
// number of individuals genotyped at both, the sums of the genotypes and of
// their squares over those individuals, and the sum of the products, are
// then popcounts of the bitwise ANDs of their planes (as PLINK does).
//
// The bytes of the FBM are decoded through its code256, so that imputed
// genotypes are used if the gen_tibble is set to use them, and are missing
// otherwise; decoded values other than 0, 1 and 2 are missing.

#include <bigstatsr/BMAcc.h>
#include <cmath>
#include <cstdint>
#include <vector>

/******************************************************************************/

// Hardware popcount: the default x86-64 target used to compile R packages does
// not include the POPCNT instruction, so the r^2 kernel is compiled twice, and
// the version with POPCNT is used on CPUs that have it (as in
// pairwise_bitpacked.h).
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#define LD_X86_POPCNT 1
#endif

#if defined(__GNUC__)
#define LD_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define LD_ALWAYS_INLINE inline
#endif

LD_ALWAYS_INLINE int64_t ld_pop(uint64_t x) {
  return static_cast<int64_t>(__builtin_popcountll(x));
}

/******************************************************************************/

// The loci packed as bit planes: for locus k, planes[k * 3 * nw + p * nw + w]
// is word w of plane p (0: genotyped, 1: genotype >= 1, 2: genotype == 2);
// and, over its genotyped individuals, their number, the sum of the genotypes
// and the sum of their squares.
struct LdPlanes {
  std::size_t nw = 0;
  std::vector<uint64_t> planes;
  std::vector<int64_t> n_valid, sum, ssq;
  const uint64_t* locus(std::size_t k) const {
    return &planes[k * 3 * nw];
  }
};

// The genotype (0, 1, 2, or 3 for missing) of each raw byte of the FBM,
// decoded through its code256
static void ld_geno_of_byte(Environment BM, uint8_t* geno_of_byte) {
  NumericVector code256 = BM["code256"];
  if (code256.size() != 256) {
    Rcpp::stop("The FBM should have a code256 of 256 values.");
  }
  for (int b = 0; b < 256; b++) {
    double g = code256[b];
    geno_of_byte[b] = (g == 0 || g == 1 || g == 2) ?
      static_cast<uint8_t>(g) : 3;
  }
}

// Pack the loci colInd (1-based columns of the FBM) of the individuals rowInd
// (1-based rows), in parallel over loci
static LdPlanes ld_pack(Environment BM, const IntegerVector& rowInd,
                        const IntegerVector& colInd, int ncores) {
  XPtr<FBM> xpBM = BM["address"];
  const unsigned char* p_mat =
    static_cast<const unsigned char*>(xpBM->matrix());
  const std::size_t nrow_fbm = xpBM->nrow();
  const std::size_t ncol_fbm = xpBM->ncol();
  const std::size_t n = rowInd.size();
  const std::size_t m = colInd.size();
  std::vector<std::size_t> rows(n), cols(m);
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
  uint8_t geno_of_byte[256];
  ld_geno_of_byte(BM, geno_of_byte);

  LdPlanes lp;
  lp.nw = (n + 63) / 64;
  lp.planes.assign(m * 3 * lp.nw, 0);
  lp.n_valid.assign(m, 0);
  lp.sum.assign(m, 0);
  lp.ssq.assign(m, 0);
  const std::size_t nw = lp.nw;
  #pragma omp parallel for num_threads(ncores) schedule(dynamic, 256)
  for (std::size_t l = 0; l < m; l++) {
    // the bytes of this locus are contiguous in the column-major FBM
    const unsigned char* p_col = p_mat + cols[l] * nrow_fbm;
    uint64_t* valid = &lp.planes[l * 3 * nw];
    uint64_t* ge1 = valid + nw;
    uint64_t* eq2 = valid + 2 * nw;
    int64_t n_valid = 0, sum = 0, ssq = 0;
    for (std::size_t i = 0; i < n; i++) {
      const uint8_t g = geno_of_byte[p_col[rows[i]]];
      if (g == 3) continue;
      const uint64_t bit = uint64_t(1) << (i & 63);
      const std::size_t w = i >> 6;
      valid[w] |= bit;
      if (g >= 1) ge1[w] |= bit;
      if (g == 2) eq2[w] |= bit;
      n_valid++;
      sum += g;
      ssq += g * g;
    }
    lp.n_valid[l] = n_valid;
    lp.sum[l] = sum;
    lp.ssq[l] = ssq;
  }
  return lp;
}

/******************************************************************************/

// The r^2 of a pair of loci, with the number of individuals genotyped at
// both. r^2 is NaN when it is undefined (no variance at either locus over
// those individuals).
struct LdR2 {
  double r2;
  int64_t n_shared;
};

// The counts over the individuals genotyped at both loci (a and b point at
// their planes); when neither locus has missing genotypes, only the sum of
// the products needs popcounts
LD_ALWAYS_INLINE LdR2 ld_r2_body(const uint64_t* a, const uint64_t* b,
                                 std::size_t nw, bool complete,
                                 int64_t n_a, int64_t s_a, int64_t q_a,
                                 int64_t s_b, int64_t q_b) {
  const uint64_t* a1 = a + nw;
  const uint64_t* a2 = a + 2 * nw;
  const uint64_t* b1 = b + nw;
  const uint64_t* b2 = b + 2 * nw;
  int64_t nj, sx, sy, sxx, syy, sxy = 0;
  if (complete) {
    for (std::size_t w = 0; w < nw; w++) {
      sxy += ld_pop(a1[w] & b1[w]) + ld_pop(a1[w] & b2[w]) +
        ld_pop(a2[w] & b1[w]) + ld_pop(a2[w] & b2[w]);
    }
    nj = n_a;
    sx = s_a;
    sxx = q_a;
    sy = s_b;
    syy = q_b;
  } else {
    nj = 0;
    int64_t sx1 = 0, sx2 = 0, sy1 = 0, sy2 = 0;
    for (std::size_t w = 0; w < nw; w++) {
      const uint64_t va = a[w], vb = b[w];
      nj += ld_pop(va & vb);
      // (the genotype planes are 0 where a genotype is missing)
      sx1 += ld_pop(a1[w] & vb);
      sx2 += ld_pop(a2[w] & vb);
      sy1 += ld_pop(b1[w] & va);
      sy2 += ld_pop(b2[w] & va);
      sxy += ld_pop(a1[w] & b1[w]) + ld_pop(a1[w] & b2[w]) +
        ld_pop(a2[w] & b1[w]) + ld_pop(a2[w] & b2[w]);
    }
    // x = (x >= 1) + (x == 2), x^2 = (x >= 1) + 3 (x == 2)
    sx = sx1 + sx2;
    sxx = sx1 + 3 * sx2;
    sy = sy1 + sy2;
    syy = sy1 + 3 * sy2;
  }
  // (exact integer covariance and variances, times nj^2)
  const double cov = static_cast<double>(nj * sxy - sx * sy);
  const double var_x = static_cast<double>(nj * sxx - sx * sx);
  const double var_y = static_cast<double>(nj * syy - sy * sy);
  LdR2 res;
  res.r2 = cov * cov / (var_x * var_y);
  res.n_shared = nj;
  return res;
}

static LdR2 ld_r2_generic(const LdPlanes& lp, std::size_t ka, std::size_t kb,
                          int64_t n) {
  const bool complete = lp.n_valid[ka] == n && lp.n_valid[kb] == n;
  return ld_r2_body(lp.locus(ka), lp.locus(kb), lp.nw, complete,
                    lp.n_valid[ka], lp.sum[ka], lp.ssq[ka], lp.sum[kb],
                    lp.ssq[kb]);
}

#ifdef LD_X86_POPCNT
// the same code, compiled to use the POPCNT instruction
__attribute__((target("popcnt")))
static LdR2 ld_r2_popcnt(const LdPlanes& lp, std::size_t ka, std::size_t kb,
                         int64_t n) {
  const bool complete = lp.n_valid[ka] == n && lp.n_valid[kb] == n;
  return ld_r2_body(lp.locus(ka), lp.locus(kb), lp.nw, complete,
                    lp.n_valid[ka], lp.sum[ka], lp.ssq[ka], lp.sum[kb],
                    lp.ssq[kb]);
}
#endif

typedef LdR2 (*LdR2Fun)(const LdPlanes&, std::size_t, std::size_t, int64_t);

// the r^2 function for this CPU
static LdR2Fun ld_r2_fun() {
#ifdef LD_X86_POPCNT
  if (__builtin_cpu_supports("popcnt")) return ld_r2_popcnt;
#endif
  return ld_r2_generic;
}

// Whether a pair of loci counts as linked: r^2 above thr; or, when fewer than
// min_shared individuals are genotyped at both loci (so that r^2 is not
// estimated), as chosen by remove_insufficient. An undefined r^2 (no
// variance over the individuals genotyped at both loci) is never above thr,
// as in PLINK and bigsnpr.
LD_ALWAYS_INLINE bool ld_linked(const LdR2& r, double thr, int64_t min_shared,
                                bool remove_insufficient) {
  if (r.n_shared < min_shared) return remove_insufficient;
  return r.r2 > thr;
}

/******************************************************************************/

// r^2 of the pairs of loci (1-based indices into colInd) given as the
// columns of `pairs` (2 x n_pairs), with the number of individuals genotyped
// at both (for tests and checks against PLINK's --r2)
// [[Rcpp::export]]
List ld_r2_pairs_cpp(Environment BM,
                     const IntegerVector& rowInd,
                     const IntegerVector& colInd,
                     const IntegerMatrix& pairs) {
  LdPlanes lp = ld_pack(BM, rowInd, colInd, 1);
  const int64_t n = rowInd.size();
  const int m = colInd.size();
  LdR2Fun r2_fun = ld_r2_fun();
  const int n_pairs = pairs.ncol();
  NumericVector r2(n_pairs);
  IntegerVector n_shared(n_pairs);
  for (int p = 0; p < n_pairs; p++) {
    const int a = pairs(0, p), b = pairs(1, p);
    if (a < 1 || a > m || b < 1 || b > m) {
      Rcpp::stop("pairs out of bounds.");
    }
    LdR2 r = r2_fun(lp, a - 1, b - 1, n);
    r2[p] = r.r2;
    n_shared[p] = static_cast<int>(r.n_shared);
  }
  return List::create(_["r2"] = r2, _["n_shared"] = n_shared);
}

/******************************************************************************/

// Clumping of the loci of one chromosome, as bigsnpr's clumping_chr()
// (src/clumping.cpp, GPL-3, Florian Privé), with r^2 over the individuals
// genotyped at both loci. The loci are visited in order of importance
// (ordInd, 1-based); a locus is kept unless it is linked (see ld_linked())
// to a more important locus that is kept, within `size` of its position
// (pos; loci of colInd in increasing position). Loci are processed in
// parallel: a thread waits for the decision on a more important locus before
// using it, so that the result is that of the sequential algorithm.
// Returns 1 for the loci kept, 0 for the others.
// [[Rcpp::export]]
IntegerVector ld_clump_cpp(Environment BM,
                           const IntegerVector& rowInd,
                           const IntegerVector& colInd,
                           const IntegerVector& ordInd,
                           const IntegerVector& rankInd,
                           const NumericVector& pos,
                           double size,
                           double thr,
                           int min_shared,
                           bool remove_insufficient,
                           int ncores) {
  const int m = colInd.size();
  if (ordInd.size() != m || rankInd.size() != m || pos.size() != m) {
    Rcpp::stop("ordInd, rankInd and pos should have one value per locus.");
  }
  LdPlanes lp = ld_pack(BM, rowInd, colInd, ncores);
  const int64_t n = rowInd.size();
  LdR2Fun r2_fun = ld_r2_fun();
  // -1: not decided yet, 0: removed, 1: kept
  std::vector<int> keep(m, -1);
  int* p_keep = keep.data();
  const int* p_rank = rankInd.begin();
  const double* p_pos = pos.begin();
  const int* p_ord = ordInd.begin();

  #pragma omp parallel num_threads(ncores)
  {
    std::vector<int> ind_to_check;
    ind_to_check.reserve(m);

    #pragma omp for schedule(dynamic, 1)
    for (int k = 0; k < m; k++) {
      const int j0 = p_ord[k] - 1;
      // the more important loci within the window that are not removed (as
      // bigsnpr's which_to_check())
      ind_to_check.clear();
      const double pos_min = p_pos[j0] - size;
      const double pos_max = p_pos[j0] + size;
      bool not_min = true, not_max = true;
      for (int d = 1; not_max || not_min; d++) {
        if (not_max) {
          const int j = j0 + d;
          not_max = (j < m) && (p_pos[j] <= pos_max);
          if (not_max && (p_rank[j0] > p_rank[j])) {
            int kj;
            #pragma omp atomic read
            kj = p_keep[j];
            if (kj != 0) ind_to_check.push_back(j);
          }
        }
        if (not_min) {
          const int j = j0 - d;
          not_min = (j >= 0) && (p_pos[j] >= pos_min);
          if (not_min && (p_rank[j0] > p_rank[j])) {
            int kj;
            #pragma omp atomic read
            kj = p_keep[j];
            if (kj != 0) ind_to_check.push_back(j);
          }
        }
      }
      const int nb_check = ind_to_check.size();
      bool keep_j0 = true;
      bool all_checked = false;
      while (!all_checked) {
        all_checked = true;
        for (int k2 = 0; k2 < nb_check; k2++) {
          const int jk = ind_to_check[k2];
          if (jk == -1) continue;
          int kj;
          #pragma omp atomic read
          kj = p_keep[jk];
          if (kj == -1) {
            // waiting for the thread deciding on locus jk
            all_checked = false;
          } else if (kj == 0) {
            // removed: no need to check it after all
            ind_to_check[k2] = -1;
          } else {
            LdR2 r = r2_fun(lp, jk, j0, n);
            if (ld_linked(r, thr, min_shared, remove_insufficient)) {
              keep_j0 = false;
              all_checked = true;
              break;
            }
            ind_to_check[k2] = -1;
          }
        }
      }
      #pragma omp atomic write
      p_keep[j0] = keep_j0;
    }
  }
  return IntegerVector(keep.begin(), keep.end());
}
