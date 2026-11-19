#ifndef PAIRWISE_BITPACKED_H
#define PAIRWISE_BITPACKED_H

/******************************************************************************/

#include <bigstatsr/BMAcc.h>
#include <vector>
#include <cstdint>
#include <algorithm>

/******************************************************************************/

// Shared kernel for the pairwise (individual x individual) statistics IBS,
// KING-robust and allele sharing (see snp_ibs.cpp, snp_king.cpp and
// snp_as.cpp), and for the counts of loci valid in both individuals of
// PLINK-style relationship matrices (gt_grm_plink.cpp).
//
// All three statistics are sums over loci of simple functions of the two
// diploid genotypes (0, 1 or 2 alternate alleles, or missing). Rather than
// recoding blocks of genotypes as dense double indicator matrices and using
// matrix products (which needs several n x block_size double matrices in
// memory and O(n^2 m) floating point operations), genotypes are packed as
// bits, 64 loci per 64 bit word, into a few "planes" per individual (e.g. one
// bit plane for "genotype is valid", one for "genotype is 2", ...). The
// contribution of 64 loci to a pair of individuals is then obtained with a
// few bitwise operations and popcounts, as done by PLINK. Counts are exact
// integers, so the results are identical to those of the matrix products.
//
// Loci are processed in chunks (up to 64 words of 64 loci at a time), small
// enough that the planes of the individuals visited by each individual (all
// of them, or those of its group) stay in cache while its pairs are visited.
// Each pair is visited once (the statistics are symmetric, or computed for
// both directions at once), and per pair counts are accumulated in the
// packed lower triangle (column-major, including the diagonal).
//
// The bytes of the FBM are decoded through its code256, so that imputed
// genotypes (bytes 4 to 6) are used if the gen_tibble is set to use them
// (code256 is then CODE_IMPUTE_PRED, see gt_set_imputed()), and treated as
// missing otherwise (CODE_012). Decoded values other than 0, 1 and 2 (NA,
// or fractional dosages) are treated as missing.

/******************************************************************************/

// Hardware popcount: the default x86-64 target used to compile R packages does
// not include the POPCNT instruction (so __builtin_popcountll would call a
// slow software routine); on x86 we therefore compile a second version of the
// inner loop for CPUs with POPCNT (virtually all of them), and choose between
// the two at run time. On other architectures (e.g. arm64) the builtin
// already maps to a fast instruction.
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#define PW_X86_POPCNT 1
#endif

#if defined(__GNUC__)
#define PW_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define PW_ALWAYS_INLINE inline
#endif

PW_ALWAYS_INLINE uint32_t pw_popcount(uint64_t x) {
  return static_cast<uint32_t>(__builtin_popcountll(x));
}

/******************************************************************************/

// Each statistic is described by a small "policy" struct with:
// - N_PLANES: the number of bit planes per individual
// - N_COUNTS: the number of counts accumulated for each pair
// - planes(g): a bitmask (bit p set if plane p is 1) for genotype g (0, 1 or
//   2); missing genotypes have no bit set (not even the valid plane)
// - add(wi, wj, c): add the counts for one word (64 loci) of individuals
//   i and j (wi, wj point at their N_PLANES words) to c

// IBS: planes are valid, (g >= 1) and (g == 2); counts are the number of loci
// valid in both individuals, and the sum of |g_i - g_j| over those loci
// (|g_i - g_j| is the number of the two thresholds g >= 1 and g == 2 at which
// the two genotypes differ).
struct PwIbs {
  static const int N_PLANES = 3;
  static const int N_COUNTS = 2;
  static uint8_t planes(int g) {
    return 1 | ((g >= 1) << 1) | ((g == 2) << 2);
  }
  static PW_ALWAYS_INLINE void add(const uint64_t* wi, const uint64_t* wj,
                                   uint32_t* c) {
    uint64_t valid = wi[0] & wj[0];
    c[0] += pw_popcount(valid);
    c[1] += pw_popcount((wi[1] ^ wj[1]) & valid) +
      pw_popcount((wi[2] ^ wj[2]) & valid);
  }
};

// IBS with PLINK 1.9's weighted correction for missing genotypes (see
// snp_ibs.cpp): the counts of PwIbs, plus the sum of the (integer) weights of
// the loci missing in both individuals. The weights of the 64 loci of the
// word are in wts (0 for the padding bits past the last locus, which are not
// valid in either individual). Loci missing in both individuals are usually
// rare, so their bits are visited one at a time, as in PLINK.
struct PwIbsWeighted {
  static const int N_PLANES = 3;
  static const int N_COUNTS = 3;
  static uint8_t planes(int g) {
    return PwIbs::planes(g);
  }
  static PW_ALWAYS_INLINE void add(const uint64_t* wi, const uint64_t* wj,
                                   uint32_t* c, const uint32_t* wts) {
    PwIbs::add(wi, wj, c);
    uint64_t missing = ~(wi[0] | wj[0]);
    while (missing) {
      c[2] += wts[__builtin_ctzll(missing)];
      missing &= missing - 1; // clear the lowest set bit
    }
  }
};

// Add the counts of word w; only PwIbsWeighted uses the weights of the loci
// of the chunk, wts (the other statistics ignore them). The offset of the
// weights of word w is only computed here, as wts is nullptr for the other
// statistics, and an offset added to a null pointer is undefined behaviour
// (even if the result is never used)
template <class Stat>
struct PwAdd {
  static PW_ALWAYS_INLINE void run(const uint64_t* wi, const uint64_t* wj,
                                   uint32_t* c, const uint32_t*,
                                   std::size_t) {
    Stat::add(wi, wj, c);
  }
};
template <>
struct PwAdd<PwIbsWeighted> {
  static PW_ALWAYS_INLINE void run(const uint64_t* wi, const uint64_t* wj,
                                   uint32_t* c, const uint32_t* wts,
                                   std::size_t w) {
    PwIbsWeighted::add(wi, wj, c, wts + w * 64);
  }
};

// KING: planes are valid, heterozygous (g == 1), homozygous reference (g == 0)
// and homozygous alternate (g == 2); counts are the number of loci at which
// both are heterozygous, the number at which they are opposite homozygotes,
// the number at which i is heterozygous and j is valid, and the number at
// which j is heterozygous and i is valid.
struct PwKing {
  static const int N_PLANES = 4;
  static const int N_COUNTS = 4;
  static uint8_t planes(int g) {
    return 1 | ((g == 1) << 1) | ((g == 0) << 2) | ((g == 2) << 3);
  }
  static PW_ALWAYS_INLINE void add(const uint64_t* wi, const uint64_t* wj,
                                   uint32_t* c) {
    c[0] += pw_popcount(wi[1] & wj[1]);
    c[1] += pw_popcount((wi[2] & wj[3]) | (wi[3] & wj[2]));
    c[2] += pw_popcount(wi[1] & wj[0]);
    c[3] += pw_popcount(wj[1] & wi[0]);
  }
};

// Allele sharing: planes are valid, homozygous reference and homozygous
// alternate; counts are the number of loci at which both are the same
// homozygote, the number at which they are opposite homozygotes (so that
// sum((g_i - 1) * (g_j - 1)) = same - opposite), and the number of loci valid
// in both.
struct PwAlleleSharing {
  static const int N_PLANES = 3;
  static const int N_COUNTS = 3;
  static uint8_t planes(int g) {
    return 1 | ((g == 0) << 1) | ((g == 2) << 2);
  }
  static PW_ALWAYS_INLINE void add(const uint64_t* wi, const uint64_t* wj,
                                   uint32_t* c) {
    c[0] += pw_popcount((wi[1] & wj[1]) | (wi[2] & wj[2]));
    c[1] += pw_popcount((wi[1] & wj[2]) | (wi[2] & wj[1]));
    c[2] += pw_popcount(wi[0] & wj[0]);
  }
};

// Valid loci: one plane (valid) and one count, the number of loci valid in
// both individuals (the denominator of PLINK 1.9's relationship matrix with
// missing genotypes, see gt_pca()).
struct PwValid {
  static const int N_PLANES = 1;
  static const int N_COUNTS = 1;
  static uint8_t planes(int) {
    return 1;
  }
  static PW_ALWAYS_INLINE void add(const uint64_t* wi, const uint64_t* wj,
                                   uint32_t* c) {
    c[0] += pw_popcount(wi[0] & wj[0]);
  }
};

/******************************************************************************/

// Counts for individual i against all individuals j with i <= j < n, over
// n_words words of the current chunk; cnt_i points at the counts of pair
// (i, i) in the packed lower triangle; wts points at the weights of the loci
// of the chunk (only used by PwIbsWeighted, nullptr otherwise).
template <class Stat>
PW_ALWAYS_INLINE void pw_row_body(std::size_t i, std::size_t n,
                                  std::size_t n_words,
                                  const uint64_t* planes, uint32_t* cnt_i,
                                  const uint32_t* wts) {
  const std::size_t stride = n_words * Stat::N_PLANES;
  const uint64_t* pi = planes + i * stride;
  for (std::size_t j = i; j < n; j++) {
    const uint64_t* pj = planes + j * stride;
    uint32_t c[Stat::N_COUNTS] = {0};
    for (std::size_t w = 0; w < n_words; w++) {
      PwAdd<Stat>::run(pi + w * Stat::N_PLANES, pj + w * Stat::N_PLANES, c,
                       wts, w);
    }
    uint32_t* cnt_ij = cnt_i + (j - i) * Stat::N_COUNTS;
    for (int k = 0; k < Stat::N_COUNTS; k++) {
      cnt_ij[k] += c[k];
    }
  }
}

// generic version (software popcount on x86 without POPCNT)
template <class Stat>
void pw_row_generic(std::size_t i, std::size_t n, std::size_t n_words,
                    const uint64_t* planes, uint32_t* cnt_i,
                    const uint32_t* wts) {
  pw_row_body<Stat>(i, n, n_words, planes, cnt_i, wts);
}

#ifdef PW_X86_POPCNT
// the same code, compiled to use the POPCNT instruction
template <class Stat>
__attribute__((target("popcnt")))
void pw_row_popcnt(std::size_t i, std::size_t n, std::size_t n_words,
                   const uint64_t* planes, uint32_t* cnt_i,
                   const uint32_t* wts) {
  pw_row_body<Stat>(i, n, n_words, planes, cnt_i, wts);
}
#endif

/******************************************************************************/

// Packing the genotypes into bit planes is a transposition: the FBM stores
// the bytes of the individuals of a locus next to each other, whereas a word
// of a plane holds 64 loci of one individual. Each byte is decoded with a
// "spread" table, which gives the bits of all planes of its genotype at once,
// plane p at bit 16 * p. Shifting the spread values of 16 consecutive loci by
// their position (0 to 15) and or-ing them gathers 16 loci of every plane in
// a single 64-bit accumulator (locus kk of plane p at bit 16 * p + kk), so
// that each genotype costs a lookup, a shift and an or, whatever the number
// of planes (at most 4).

// fully unroll the next loop (gcc and clang; other compilers ignore it)
#if defined(__GNUC__)
#define PW_UNROLL_16 _Pragma("GCC unroll 16")
#else
#define PW_UNROLL_16
#endif

// Pack one word (64 loci) of one individual: col_ptr[k] points at the column
// (in the FBM) of locus k of the word, r is the row of the individual in the
// FBM, spread is the table above, and loci_mask clears the bits of the
// loci past the end of the last word. The N_PLANES words are written to
// p_word.
template <class Stat>
PW_ALWAYS_INLINE void pw_pack_word(const unsigned char* const* col_ptr,
                                   std::size_t r, const uint64_t* spread,
                                   uint64_t loci_mask, uint64_t* p_word) {
  uint64_t acc[4];
  for (int q = 0; q < 4; q++) {
    const unsigned char* const* cp = col_ptr + 16 * q;
    uint64_t a = 0;
    PW_UNROLL_16
    for (int kk = 0; kk < 16; kk++) {
      a |= spread[cp[kk][r]] << kk;
    }
    acc[q] = a;
  }
  // the 16 bits of plane p of each accumulator are the q-th quarter of the
  // word of that plane
  for (int p = 0; p < Stat::N_PLANES; p++) {
    uint64_t word = 0;
    for (int q = 0; q < 4; q++) {
      word |= ((acc[q] >> (16 * p)) & 0xFFFFull) << (16 * q);
    }
    p_word[p] = word & loci_mask;
  }
}

/******************************************************************************/

// Accumulate the counts of statistic Stat for all pairs of the selected
// individuals (rowInd) over the selected loci (colInd).
// Returns the counts for the packed lower triangle (pairs (j, i) with j >= i,
// column-major), N_COUNTS consecutive values per pair; tri_start[i] is the
// index of pair (i, i).
// If row_end is not empty, only the pairs within groups of consecutive
// individuals are counted: row_end[i] is the (0-based, exclusive) end of the
// group of individual i in rowInd, and individual i is only paired with
// individuals j with i <= j < row_end[i] (so the packed triangles of the
// groups follow one another).
// locus_wts (only for PwIbsWeighted) holds one integer weight per selected
// locus (colInd).
template <class Stat>
std::vector<uint32_t> pw_bitpacked_counts(Environment BM,
                                          const IntegerVector& rowInd,
                                          const IntegerVector& colInd,
                                          std::vector<std::size_t>& tri_start,
                                          int ncores,
                                          const std::vector<std::size_t>&
                                            row_end =
                                            std::vector<std::size_t>(),
                                          const std::vector<uint32_t>&
                                            locus_wts =
                                            std::vector<uint32_t>()) {

  XPtr<FBM> xpBM = BM["address"];
  // raw (undecoded) bytes of the FBM, stored column-major
  const unsigned char* p_mat =
    static_cast<const unsigned char*>(xpBM->matrix());
  std::size_t nrow_fbm = xpBM->nrow();
  std::size_t ncol_fbm = xpBM->ncol();

  std::size_t n = rowInd.size(); // number of individuals
  std::size_t m = colInd.size(); // number of loci

  // 0-based offsets of rows and columns, checking the indices
  std::vector<std::size_t> rows(n), cols(m);
  for (std::size_t i = 0; i < n; i++) {
    if (rowInd[i] < 1 || static_cast<std::size_t>(rowInd[i]) > nrow_fbm) {
      Rcpp::stop("Subscript out of bounds.");
    }
    rows[i] = rowInd[i] - 1;
  }
  for (std::size_t j = 0; j < m; j++) {
    if (colInd[j] < 1 || static_cast<std::size_t>(colInd[j]) > ncol_fbm) {
      Rcpp::stop("Subscript out of bounds.");
    }
    cols[j] = colInd[j] - 1;
  }

  // the end of the pairs of each individual (all individuals, or the end of
  // its group)
  std::vector<std::size_t> pair_end(n, n);
  if (!row_end.empty()) {
    if (row_end.size() != n) {
      Rcpp::stop("row_end should have one value per individual.");
    }
    for (std::size_t i = 0; i < n; i++) {
      if (row_end[i] <= i || row_end[i] > n) {
        Rcpp::stop("row_end should be after each individual.");
      }
      pair_end[i] = row_end[i];
    }
  }

  // start of each column of the packed lower triangle
  tri_start.resize(n);
  std::size_t n_pairs = 0;
  for (std::size_t i = 0; i < n; i++) {
    tri_start[i] = n_pairs;
    n_pairs += pair_end[i] - i;
  }
  std::vector<uint32_t> counts(n_pairs * Stat::N_COUNTS, 0);

  // weights of the loci, padded with zeros to a whole number of words (so
  // that the padding bits of the last word, which look missing in both
  // individuals, have no weight); chunks start at multiples of 64 loci, so
  // the weights of a chunk start at p_wts + l0
  std::vector<uint32_t> wts_padded;
  if (!locus_wts.empty()) {
    if (locus_wts.size() != m) {
      Rcpp::stop("locus_wts should have one value per locus.");
    }
    wts_padded.assign((m + 63) / 64 * 64, 0);
    std::copy(locus_wts.begin(), locus_wts.end(), wts_padded.begin());
  }
  const uint32_t* p_wts = wts_padded.empty() ? nullptr : wts_padded.data();

  // planes for each possible raw byte, decoded through the code256 of the
  // FBM (NA, and values other than 0, 1 and 2, are missing)
  NumericVector code256 = BM["code256"];
  if (code256.size() != 256) {
    Rcpp::stop("The FBM should have a code256 of 256 values.");
  }
  // (spread table: plane p of the genotype at bit 16 * p, see pw_pack_word())
  static_assert(Stat::N_PLANES <= 4, "At most 4 planes can be packed.");
  uint64_t spread[256];
  for (int b = 0; b < 256; b++) {
    double g = code256[b];
    uint64_t bits = (g == 0 || g == 1 || g == 2) ?
      Stat::planes(static_cast<int>(g)) : 0;
    spread[b] = 0;
    for (int p = 0; p < Stat::N_PLANES; p++) {
      spread[b] |= ((bits >> p) & 1) << (16 * p);
    }
  }

  // individuals are packed in the order of their rows in the FBM, so that
  // the bytes read for a block of individuals are close together in each
  // column (rowInd is in the order of the groups, not of the FBM); their
  // planes are stored at their position in rowInd
  std::vector<std::size_t> pack_order(n);
  for (std::size_t i = 0; i < n; i++) {
    pack_order[i] = i;
  }
  std::sort(pack_order.begin(), pack_order.end(),
            [&rows](std::size_t a, std::size_t b) {
              return rows[a] < rows[b];
            });

  // number of 64-bit words per chunk: keep the planes of the individuals
  // visited by each individual (all of them, or those of its group) within
  // ~256KB, so that they stay in cache whilst visiting the pairs, and the
  // planes of all individuals within ~64MB; but use at least 4 and at most
  // 64 words per chunk
  std::size_t max_span = 1;
  for (std::size_t i = 0; i < n; i++) {
    max_span = std::max(max_span, pair_end[i] - i);
  }
  std::size_t chunk_words = (256 * 1024) / (8 * Stat::N_PLANES * max_span);
  chunk_words = std::min(chunk_words, (64 * 1024 * 1024) /
    (8 * Stat::N_PLANES * std::max<std::size_t>(n, 1)));
  chunk_words =
    std::min<std::size_t>(64, std::max<std::size_t>(4, chunk_words));
  const std::size_t chunk_loci = chunk_words * 64;
  // planes, for each individual: chunk_words groups of N_PLANES words
  std::vector<uint64_t> planes(n * chunk_words * Stat::N_PLANES);
  uint64_t* p_planes = planes.data();
  uint32_t* p_counts = counts.data();

  // use the hardware popcount when available
  bool use_popcnt = false;
#ifdef PW_X86_POPCNT
  use_popcnt = __builtin_cpu_supports("popcnt");
#endif

  #pragma omp parallel num_threads(ncores)
  {
    for (std::size_t l0 = 0; l0 < m; l0 += chunk_loci) {
      std::size_t n_loci = std::min(chunk_loci, m - l0);
      std::size_t n_words = (n_loci + 63) / 64;
      std::size_t stride = n_words * Stat::N_PLANES;

      // pack the genotypes of this chunk into bit planes, by words of 64
      // loci and blocks of 64 individuals (in the order of the FBM); each
      // word and block is a unit of work, and the blocks of a word follow
      // one another, so that the 64 columns of a word are read from start
      // to end (rather than jumping to other columns after 64 individuals)
      std::size_t n_blocks = (n + 63) / 64;
      #pragma omp for schedule(static)
      for (std::size_t unit = 0; unit < n_blocks * n_words; unit++) {
        std::size_t i0 = (unit % n_blocks) * 64;
        std::size_t w = unit / n_blocks;
        std::size_t n_ind = std::min<std::size_t>(64, n - i0);
        std::size_t n_bits = std::min<std::size_t>(64, n_loci - w * 64);
        // the columns of the loci of this word (the loci past the end of the
        // last word point at its first column, and are masked out)
        const unsigned char* col_ptr[64];
        for (std::size_t k = 0; k < 64; k++) {
          std::size_t l = l0 + w * 64 + (k < n_bits ? k : 0);
          col_ptr[k] = p_mat + cols[l] * nrow_fbm;
        }
        uint64_t loci_mask =
          (n_bits == 64) ? ~0ull : ((1ull << n_bits) - 1);
        for (std::size_t ii = 0; ii < n_ind; ii++) {
          std::size_t idx = pack_order[i0 + ii];
          pw_pack_word<Stat>(col_ptr, rows[idx], spread, loci_mask,
                             p_planes + idx * stride + w * Stat::N_PLANES);
        }
      }
      // (implicit barrier: all planes are ready)

      // visit all pairs (j >= i); rows of the triangle have different
      // lengths, so use a dynamic schedule
      #pragma omp for schedule(dynamic, 8)
      for (std::size_t i = 0; i < n; i++) {
        uint32_t* cnt_i = p_counts + tri_start[i] * Stat::N_COUNTS;
        std::size_t end_i = pair_end[i];
        const uint32_t* wts_chunk = p_wts ? p_wts + l0 : nullptr;
#ifdef PW_X86_POPCNT
        if (use_popcnt) {
          pw_row_popcnt<Stat>(i, end_i, n_words, p_planes, cnt_i, wts_chunk);
        } else {
          pw_row_generic<Stat>(i, end_i, n_words, p_planes, cnt_i, wts_chunk);
        }
#else
        pw_row_generic<Stat>(i, end_i, n_words, p_planes, cnt_i, wts_chunk);
#endif
      }
      // (implicit barrier: the planes can be overwritten by the next chunk)
    }
  }

  return counts;
}

// get a pointer to the data of a double FBM (as used for the accumulators)
inline double* pw_fbm_double_ptr(Environment BM, std::size_t n) {
  XPtr<FBM_RW> xpBM = BM["address_rw"];
  if (xpBM->matrix_type() != 8) {
    Rcpp::stop("Accumulators should be 'double' FBMs.");
  }
  if (xpBM->nrow() != n || xpBM->ncol() != n) {
    Rcpp::stop("Accumulators should be n x n FBMs.");
  }
  return static_cast<double*>(xpBM->matrix());
}

/******************************************************************************/

#endif // PAIRWISE_BITPACKED_H
