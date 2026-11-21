#ifndef PAIRWISE_BITPACKED_H
#define PAIRWISE_BITPACKED_H

/******************************************************************************/

#include <bigstatsr/BMAcc.h>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cstdlib>

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

// AVX2 (x86-64 only): KING counts for 4 words at a time, with a vectorised
// popcount (a lookup of each half byte with vpshufb, Mula, Kurz and Lemire
// 2018, https://doi.org/10.1093/comjnl/bxx046), compiled with the target
// attribute and chosen at run time on CPUs with AVX2 (no compiler flags
// needed)
#if defined(__GNUC__) && defined(__x86_64__)
#define PW_X86_AVX2 1
#include <immintrin.h>
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

// Layout of the planes of a chunk, for each individual: with blk = 1 (the
// scalar kernels), the N_PLANES words of word 0, then those of word 1, etc.;
// with blk = 4 (the AVX2 kernel), blocks of 4 words (256 loci), and within a
// block the 4 words of plane 0, then the 4 words of plane 1, etc. (so that
// SIMD code can load the 4 words of a plane at once). Word w of plane p is at
// (w / blk) * blk * N_PLANES + blk * p + w % blk. The number of words of an
// individual, ps, is the number of words of the chunk rounded up to a
// multiple of blk, and the padding words are 0. (The scalar kernels read the
// blocked layout of blk = 4 at half the speed, hence the two layouts.)
template <class Stat>
PW_ALWAYS_INLINE std::size_t pw_word_offset(std::size_t w, int p,
                                            std::size_t blk) {
  return (w / blk) * blk * Stat::N_PLANES + blk * p + w % blk;
}

// Counts for individual i against the individuals j with j0 <= j < j1 (and
// j >= i), over the n_words words of the current chunk (layout blk = 1, see
// above); cnt_i points at the counts of pair (i, i) in the packed lower
// triangle; wts points at the weights of the loci of the chunk (only used by
// PwIbsWeighted, nullptr otherwise).
template <class Stat>
PW_ALWAYS_INLINE void pw_row_body(std::size_t i, std::size_t j0,
                                  std::size_t j1, std::size_t n_words,
                                  const uint64_t* planes, uint32_t* cnt_i,
                                  const uint32_t* wts) {
  const std::size_t stride = n_words * Stat::N_PLANES;
  const uint64_t* pi = planes + i * stride;
  for (std::size_t j = j0; j < j1; j++) {
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
void pw_row_generic(std::size_t i, std::size_t j0, std::size_t j1,
                    std::size_t n_words, const uint64_t* planes,
                    uint32_t* cnt_i, const uint32_t* wts) {
  pw_row_body<Stat>(i, j0, j1, n_words, planes, cnt_i, wts);
}

#ifdef PW_X86_POPCNT
// the same code, compiled to use the POPCNT instruction
template <class Stat>
__attribute__((target("popcnt")))
void pw_row_popcnt(std::size_t i, std::size_t j0, std::size_t j1,
                   std::size_t n_words, const uint64_t* planes,
                   uint32_t* cnt_i, const uint32_t* wts) {
  pw_row_body<Stat>(i, j0, j1, n_words, planes, cnt_i, wts);
}
#endif

// SIMD versions of the row kernel, for the statistics that have one (only
// KING for now): PwAvx2<Stat>::available says whether there is one
template <class Stat>
struct PwAvx2 {
  static const bool available = false;
  static void run(std::size_t, std::size_t, std::size_t, std::size_t,
                  const uint64_t*, uint32_t*) {}
};

#ifdef PW_X86_AVX2
// the number of set bits of each byte of v (a lookup of each half byte)
#define PW_POPCNT8_AVX2(v)                                                \
  _mm256_add_epi8(                                                        \
    _mm256_shuffle_epi8(lut, _mm256_and_si256((v), low)),                 \
    _mm256_shuffle_epi8(lut,                                              \
                        _mm256_and_si256(_mm256_srli_epi16((v), 4), low)))

// sum of the 4 64-bit lanes of the byte sums of v
#define PW_HSUM8_AVX2(v, out)                                             \
  do {                                                                    \
    __m256i s_ = _mm256_sad_epu8((v), zero);                              \
    __m128i t_ = _mm_add_epi64(_mm256_castsi256_si128(s_),                \
                               _mm256_extracti128_si256(s_, 1));          \
    (out) = static_cast<uint32_t>(_mm_cvtsi128_si64(t_) +                 \
                                  _mm_extract_epi64(t_, 1));              \
  } while (0)

// KING counts (see PwKing) of individual i against individuals j0 <= j < j1,
// 4 words at a time; the byte counters cannot overflow, as a chunk has at
// most 64 words (16 blocks, at most 8 per byte and block)
__attribute__((target("avx2")))
inline void pw_row_king_avx2(std::size_t i, std::size_t j0, std::size_t j1,
                             std::size_t ps, const uint64_t* planes,
                             uint32_t* cnt_i) {
  const std::size_t stride = ps * PwKing::N_PLANES;
  const std::size_t n_blocks = ps / 4;
  const uint64_t* pi = planes + i * stride;
  const __m256i lut = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3,
                                       2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3,
                                       1, 2, 2, 3, 2, 3, 3, 4);
  const __m256i low = _mm256_set1_epi8(0x0f);
  const __m256i zero = _mm256_setzero_si256();
  for (std::size_t j = j0; j < j1; j++) {
    const uint64_t* pj = planes + j * stride;
    __m256i c0 = zero, c1 = zero, c2 = zero, c3 = zero;
    for (std::size_t b = 0; b < n_blocks; b++) {
      const __m256i* bi = reinterpret_cast<const __m256i*>(pi + b * 16);
      const __m256i* bj = reinterpret_cast<const __m256i*>(pj + b * 16);
      // planes: valid, heterozygous, homozygous reference and alternate
      const __m256i vi = _mm256_loadu_si256(bi);
      const __m256i hi = _mm256_loadu_si256(bi + 1);
      const __m256i ri = _mm256_loadu_si256(bi + 2);
      const __m256i ai = _mm256_loadu_si256(bi + 3);
      const __m256i vj = _mm256_loadu_si256(bj);
      const __m256i hj = _mm256_loadu_si256(bj + 1);
      const __m256i rj = _mm256_loadu_si256(bj + 2);
      const __m256i aj = _mm256_loadu_si256(bj + 3);
      const __m256i hethet = _mm256_and_si256(hi, hj);
      const __m256i opposite = _mm256_or_si256(_mm256_and_si256(ri, aj),
                                               _mm256_and_si256(ai, rj));
      const __m256i het_i = _mm256_and_si256(hi, vj);
      const __m256i het_j = _mm256_and_si256(hj, vi);
      c0 = _mm256_add_epi8(c0, PW_POPCNT8_AVX2(hethet));
      c1 = _mm256_add_epi8(c1, PW_POPCNT8_AVX2(opposite));
      c2 = _mm256_add_epi8(c2, PW_POPCNT8_AVX2(het_i));
      c3 = _mm256_add_epi8(c3, PW_POPCNT8_AVX2(het_j));
    }
    uint32_t* cnt_ij = cnt_i + (j - i) * PwKing::N_COUNTS;
    uint32_t s;
    PW_HSUM8_AVX2(c0, s);
    cnt_ij[0] += s;
    PW_HSUM8_AVX2(c1, s);
    cnt_ij[1] += s;
    PW_HSUM8_AVX2(c2, s);
    cnt_ij[2] += s;
    PW_HSUM8_AVX2(c3, s);
    cnt_ij[3] += s;
  }
}

template <>
struct PwAvx2<PwKing> {
  static const bool available = true;
  static void run(std::size_t i, std::size_t j0, std::size_t j1,
                  std::size_t ps, const uint64_t* planes, uint32_t* cnt_i) {
    pw_row_king_avx2(i, j0, j1, ps, planes, cnt_i);
  }
};
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
// p_word[0], p_word[blk], ... (see pw_word_offset()).
template <class Stat>
PW_ALWAYS_INLINE void pw_pack_word(const unsigned char* const* col_ptr,
                                   std::size_t r, const uint64_t* spread,
                                   uint64_t loci_mask, std::size_t blk,
                                   uint64_t* p_word) {
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
  PW_UNROLL_16
  for (int p = 0; p < Stat::N_PLANES; p++) {
    uint64_t word = 0;
    PW_UNROLL_16
    for (int q = 0; q < 4; q++) {
      word |= ((acc[q] >> (16 * p)) & 0xFFFFull) << (16 * q);
    }
    p_word[blk * p] = word & loci_mask;
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

  // number of 64-bit words per chunk: 64 (4096 loci), unless the planes of
  // all individuals would exceed ~64MB; but at least 4 words. The pairs are
  // visited by tiles (below), so that the planes of the individuals visited
  // stay in cache whatever the size of the chunk; long chunks amortise the
  // update of the counts of each pair (which, for many individuals, do not
  // fit in cache)
  std::size_t chunk_words = (64 * 1024 * 1024) /
    (8 * Stat::N_PLANES * std::max<std::size_t>(n, 1));
  chunk_words =
    std::min<std::size_t>(64, std::max<std::size_t>(4, chunk_words));
  const std::size_t chunk_loci = chunk_words * 64;
  // planes, for each individual: chunk_words (rounded up to a multiple of 4
  // for the AVX2 layout) words of each of the N_PLANES planes (see
  // pw_word_offset())
  const std::size_t chunk_words_pad = (chunk_words + 3) / 4 * 4;
  std::vector<uint64_t> planes(n * chunk_words_pad * Stat::N_PLANES);
  uint64_t* p_planes = planes.data();
  uint32_t* p_counts = counts.data();

  // the pairs are visited by tiles of individuals: a unit of work is the
  // pairs (i, j) with i in tile ti and j in tile tj (ti <= tj, j >= i), so
  // that the planes of the individuals of tile tj (~128KB) stay in cache
  // whilst the individuals of tile ti are compared with them; units
  // without pairs (e.g. across groups) are skipped
  const std::size_t tile = std::min<std::size_t>(256, std::max<std::size_t>(
    16, (128 * 1024) / (8 * chunk_words * Stat::N_PLANES)));
  const std::size_t n_tiles = (n + tile - 1) / tile;
  std::vector<std::size_t> tile_max_end(n_tiles, 0);
  for (std::size_t i = 0; i < n; i++) {
    tile_max_end[i / tile] = std::max(tile_max_end[i / tile], pair_end[i]);
  }
  std::vector<std::size_t> unit_ti, unit_tj;
  for (std::size_t tj = 0; tj < n_tiles; tj++) {
    for (std::size_t ti = 0; ti <= tj; ti++) {
      if (tile_max_end[ti] > tj * tile) {
        unit_ti.push_back(ti);
        unit_tj.push_back(tj);
      }
    }
  }
  const std::size_t n_units = unit_ti.size();

  // use the hardware popcount when available
  bool use_popcnt = false;
#ifdef PW_X86_POPCNT
  use_popcnt = __builtin_cpu_supports("popcnt");
#endif
  // and the AVX2 kernel (if the statistic has one); the environment
  // variable TIDYPOPGEN_NO_AVX2 turns it off (to test the scalar kernel on
  // CPUs with AVX2, and for timings)
  bool use_avx2 = false;
#ifdef PW_X86_AVX2
  use_avx2 = PwAvx2<Stat>::available && __builtin_cpu_supports("avx2") &&
    std::getenv("TIDYPOPGEN_NO_AVX2") == nullptr;
#endif
  // the layout of the planes (see pw_word_offset())
  const std::size_t blk = use_avx2 ? 4 : 1;

  #pragma omp parallel num_threads(ncores)
  {
    for (std::size_t l0 = 0; l0 < m; l0 += chunk_loci) {
      std::size_t n_loci = std::min(chunk_loci, m - l0);
      std::size_t n_words = (n_loci + 63) / 64;
      // words per plane (a multiple of blk) and per individual
      std::size_t ps = (n_words + blk - 1) / blk * blk;
      std::size_t stride = ps * Stat::N_PLANES;

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
        // prefetch the bytes of the next block of individuals (the next
        // unit of this thread) in these columns: the bytes of a block are
        // fetched from memory column by column, a cache line from each of 64
        // columns, which the hardware does not prefetch well; the first and
        // last rows of the block cover it if its rows are consecutive. They
        // are prefetched into L2 (locality 1), which allows more misses in
        // flight than L1 (prefetching into L1, or two blocks ahead, was
        // slower)
#if defined(__GNUC__)
        if (i0 + 64 < n) {
          std::size_t r_first = rows[pack_order[i0 + 64]];
          std::size_t r_last =
            rows[pack_order[std::min<std::size_t>(i0 + 127, n - 1)]];
          for (std::size_t k = 0; k < n_bits; k++) {
            __builtin_prefetch(col_ptr[k] + r_first, 0, 1);
            __builtin_prefetch(col_ptr[k] + r_last, 0, 1);
          }
        }
#endif
        // pack each individual of the block
        for (std::size_t ii = 0; ii < n_ind; ii++) {
          std::size_t idx = pack_order[i0 + ii];
          pw_pack_word<Stat>(col_ptr, rows[idx], spread, loci_mask, blk,
                             p_planes + idx * stride +
                               pw_word_offset<Stat>(w, 0, blk));
        }
      }
      // (implicit barrier) the padding words of the last block are 0
      if (ps > n_words) {
        #pragma omp for schedule(static)
        for (std::size_t idx = 0; idx < n; idx++) {
          for (std::size_t w = n_words; w < ps; w++) {
            for (int p = 0; p < Stat::N_PLANES; p++) {
              p_planes[idx * stride + pw_word_offset<Stat>(w, p, blk)] = 0;
            }
          }
        }
      }
      // (implicit barrier: all planes are ready)

      // visit all pairs (j >= i), tile by tile; tiles on the diagonal have
      // half the pairs, so use a dynamic schedule
      const uint32_t* wts_chunk = p_wts ? p_wts + l0 : nullptr;
      #pragma omp for schedule(dynamic, 1)
      for (std::size_t u = 0; u < n_units; u++) {
        std::size_t i0 = unit_ti[u] * tile;
        std::size_t i1 = std::min(n, i0 + tile);
        std::size_t j0 = unit_tj[u] * tile;
        std::size_t j1 = std::min(n, j0 + tile);
        for (std::size_t i = i0; i < i1; i++) {
          std::size_t start = std::max(i, j0);
          std::size_t end = std::min(pair_end[i], j1);
          if (start >= end) {
            continue;
          }
          uint32_t* cnt_i = p_counts + tri_start[i] * Stat::N_COUNTS;
          if (use_avx2) {
            PwAvx2<Stat>::run(i, start, end, ps, p_planes, cnt_i);
            continue;
          }
#ifdef PW_X86_POPCNT
          if (use_popcnt) {
            pw_row_popcnt<Stat>(i, start, end, n_words, p_planes, cnt_i,
                                wts_chunk);
          } else {
            pw_row_generic<Stat>(i, start, end, n_words, p_planes, cnt_i,
                                 wts_chunk);
          }
#else
          pw_row_generic<Stat>(i, start, end, n_words, p_planes, cnt_i,
                               wts_chunk);
#endif
        }
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
