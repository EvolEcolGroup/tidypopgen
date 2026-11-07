#ifndef GT_IND_COUNTS_H
#define GT_IND_COUNTS_H

/******************************************************************************/

#include "fbm_byte_acc.h"
#include <vector>
#include <cstdint>
#include <algorithm>

/******************************************************************************/

// Shared kernel for per-individual (row-wise) counts of heterozygous and
// missing genotypes, used by gt_ind_hetero() and gt_ind_hetero_polyploid().
//
// Rather than decoding every genotype through the code256 accessor (a double
// lookup plus data-dependent branches, which the compiler cannot keep in
// registers because of the accessor's indirections), each raw byte of the FBM
// is mapped through a small lookup table (built once from code256, and from
// each individual's ploidy) to a packed 64 bit increment: the number of
// heterozygotes in the low 32 bits and the number of missing genotypes in the
// high 32 bits. The inner loop is then a single load and add per genotype,
// with no branches. Individuals that need different tables (e.g. different
// ploidy) are given different lut "classes", each a block of 256 entries.
// The genotypes are read through FBMByteAcc (see fbm_byte_acc.h), which
// resolves each locus to a pointer to its column without copying the column
// indices.

// Pack a heterozygote and a missing indicator into a single increment
inline uint64_t gt_pack_het_na(bool het, bool na) {
  return static_cast<uint64_t>(het) | (static_cast<uint64_t>(na) << 32);
}

// The inner loop over individuals, for a group of N_COLS loci (p_cols points at
// the start of each locus in the FBM). Processing several loci at once means
// that each individual's counter is read and written once per group rather
// than once per genotype. The two template flags drop the extra loads that are
// not needed in the most common case (all rows of the FBM in their original
// order, and a single lookup table, as for diploid data).
template <int N_COLS, bool all_rows, bool single_class>
inline void gt_ind_het_na_cols(const unsigned char* const* p_cols,
                               const uint64_t* p_lut,
                               const std::size_t* p_rows,
                               const std::size_t* p_class,
                               uint64_t* p_counts,
                               std::size_t n) {
  for (std::size_t i = 0; i < n; i++) {
    std::size_t row = all_rows ? i : p_rows[i];
    const uint64_t* this_lut = single_class ? p_lut : p_lut + p_class[i];
    uint64_t incr = 0;
    for (int k = 0; k < N_COLS; k++) {
      incr += this_lut[p_cols[k][row]];
    }
    p_counts[i] += incr;
  }
}

// dispatch on the template flags (for N_COLS loci)
template <int N_COLS>
inline void gt_ind_het_na_dispatch(bool all_rows, bool single_class,
                                   const unsigned char* const* p_cols,
                                   const uint64_t* p_lut,
                                   const std::size_t* p_rows,
                                   const std::size_t* p_class,
                                   uint64_t* p_counts,
                                   std::size_t n) {
  if (all_rows && single_class) {
    gt_ind_het_na_cols<N_COLS, true, true>(p_cols, p_lut, p_rows, p_class,
                                           p_counts, n);
  } else if (all_rows) {
    gt_ind_het_na_cols<N_COLS, true, false>(p_cols, p_lut, p_rows, p_class,
                                            p_counts, n);
  } else if (single_class) {
    gt_ind_het_na_cols<N_COLS, false, true>(p_cols, p_lut, p_rows, p_class,
                                            p_counts, n);
  } else {
    gt_ind_het_na_cols<N_COLS, false, false>(p_cols, p_lut, p_rows, p_class,
                                             p_counts, n);
  }
}

// Count heterozygotes and missing genotypes for each individual.
// - lut: n_classes * 256 packed increments (see gt_pack_het_na())
// - row_class: for each selected individual, the offset (class * 256) of its
//   block in lut (all 0 when there is a single class)
// Returns a 2 x n integer matrix (row 1: n_het, row 2: n_na).
inline IntegerMatrix gt_ind_het_na_counts(Environment BM,
                                          const IntegerVector& rowInd,
                                          const IntegerVector& colInd,
                                          const std::vector<uint64_t>& lut,
                                          const std::vector<std::size_t>& row_class,
                                          int ncores) {

  // raw (undecoded) bytes of the selected rows and columns of the FBM
  // (indices are checked against the dimensions of the FBM)
  FBMByteAcc acc(BM, rowInd, colInd);
  const std::vector<std::size_t>& rows = acc.rows(); // 0-based

  std::size_t n = acc.nrow(); // number of individuals
  std::size_t m = acc.ncol(); // number of loci

  // are we using all the rows of the FBM, in their original order?
  bool all_rows = (n == acc.nrow_fbm());
  for (std::size_t i = 0; all_rows && i < n; i++) {
    all_rows = (rows[i] == i);
  }
  // do all individuals share the same lookup table?
  bool single_class = (lut.size() == 256);

  // packed counts per individual (het in the low, na in the high 32 bits)
  std::vector<uint64_t> counts(n, 0);

  // the outer loop is parallelised over loci (columns are contiguous in the
  // column-major FBM), but the output is indexed by individual, so each thread
  // accumulates into its own counters, which are merged once per thread
  #pragma omp parallel num_threads(ncores)
  {
    std::vector<uint64_t> counts_local(n, 0);
    uint64_t* p_counts = counts_local.data();
    const uint64_t* p_lut = lut.data();
    const std::size_t* p_rows = rows.data();
    const std::size_t* p_class = row_class.data();

    // loci are processed in groups of 4 (the last group can be smaller)
    const std::size_t n_groups = (m + 3) / 4;
    #pragma omp for schedule(static)
    for (std::size_t g = 0; g < n_groups; g++) {
      std::size_t j0 = g * 4;
      std::size_t n_cols = std::min<std::size_t>(4, m - j0);
      // start of each locus of this group in the FBM
      const unsigned char* p_cols[4];
      for (std::size_t k = 0; k < n_cols; k++) {
        p_cols[k] = acc.col(j0 + k);
      }
      if (n_cols == 4) {
        gt_ind_het_na_dispatch<4>(all_rows, single_class, p_cols, p_lut,
                                  p_rows, p_class, p_counts, n);
      } else {
        // the remaining (up to 3) loci, one at a time
        for (std::size_t k = 0; k < n_cols; k++) {
          gt_ind_het_na_dispatch<1>(all_rows, single_class, p_cols + k, p_lut,
                                    p_rows, p_class, p_counts, n);
        }
      }
    }

    #pragma omp critical
    {
      for (std::size_t i = 0; i < n; i++) {
        counts[i] += counts_local[i];
      }
    }
  }

  // unpack the counts into a matrix of 2 rows, n_het and n_na
  IntegerMatrix het_counts(2, n);
  for (std::size_t i = 0; i < n; i++) {
    het_counts(0, i) = static_cast<int>(counts[i] & 0xFFFFFFFFu);
    het_counts(1, i) = static_cast<int>(counts[i] >> 32);
  }
  return het_counts;
}

/******************************************************************************/

#endif // GT_IND_COUNTS_H
