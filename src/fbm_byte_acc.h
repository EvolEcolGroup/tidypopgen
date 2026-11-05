#ifndef TIDYPOPGEN_FBM_BYTE_ACC_H
#define TIDYPOPGEN_FBM_BYTE_ACC_H

/******************************************************************************/
// Fast, read-only access to the genotypes of an FBM.code256, used by the
// per-locus kernels (allele frequencies, missingness, pi, HWE, ...).
//
// bigstatsr's SubBMCode256Acc decodes every genotype through two index
// vectors (rowInd and colInd, both copied to std::vector<size_t>, i.e. 8 bytes
// per locus) and a bounds-checked Rcpp NumericVector lookup, returning a
// double; kernels then accumulate those doubles into Rcpp matrices inside the
// OpenMP loop. On HGDP (1043 x 643733) this runs at ~110 M genotypes/s, about
// 15x slower than simply reading the bytes. Here instead:
// - each locus is resolved once to a pointer to the start of its column, and
//   the column indices are not copied;
// - individuals are reordered by group (stable, so the order within a group is
//   kept) and split into runs that share a group and a lookup table, and runs
//   of consecutive rows are read directly from the column pointer;
// - code256 is turned into integer lookup tables on the raw bytes, built at
//   call time (so imputed codes keep working), where each entry packs several
//   counters (e.g. alt allele count and valid allele count) into one 64-bit
//   word, so each genotype costs a single integer addition;
// - no Rcpp object is touched inside the OpenMP regions: totals are kept in
//   thread-local buffers and results are written through raw pointers.
//
// Kernels either loop over all loci with lut_sums() (push interface, with a
// callback that stores the results of each locus), or call
// LutSums::locus() from their own loops (pull interface). alt_valid_counts()
// wraps the common case of alternate and valid allele counts by group.
/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

/******************************************************************************/

// Accessor to the raw bytes of a subset (rows and columns) of an FBM.code256
class FBMByteAcc {
public:
  FBMByteAcc(Environment BM,
             const IntegerVector& rowInd,
             const IntegerVector& colInd) {
    XPtr<FBM> xpBM = BM["address"];
    _pMat = static_cast<const unsigned char*>(xpBM->matrix());
    _nrow_fbm = xpBM->nrow();
    // 0-based row indices (checking bounds)
    _rows = vec_int_to_size(rowInd, _nrow_fbm, 1);
    // check the bounds of the (1-based) column indices, but keep using the
    // integer vector itself rather than a copy
    size_t ncol_fbm = xpBM->ncol();
    _ncol = colInd.size();
    _colInd = colInd.begin();
    for (size_t j = 0; j < _ncol; j++) {
      myassert_bounds(static_cast<size_t>(_colInd[j] - 1), ncol_fbm);
    }
    // copy the code used to decode each byte into a plain array
    NumericVector code256 = BM["code256"];
    for (int b = 0; b < 256; b++) _code256[b] = code256[b];
  }

  size_t nrow() const { return _rows.size(); }
  size_t ncol() const { return _ncol; }
  const std::vector<size_t>& rows() const { return _rows; }
  // value decoded from byte b (NA for missing genotypes)
  double code(int b) const { return _code256[b]; }

  // pointer to the first byte (first row of the FBM) of the j-th locus
  inline const unsigned char* col(size_t j) const {
    return _pMat + _nrow_fbm * static_cast<size_t>(_colInd[j] - 1);
  }

  // TRUE if every non-missing value of code256 is an integer between 0 and
  // max_value, so that lookup tables of integer counts give exact results
  bool has_int_codes(double max_value) const {
    for (int b = 0; b < 256; b++) {
      double x = _code256[b];
      if (!ISNAN(x) && !(x >= 0 && x <= max_value && x == (int) x)) {
        return false;
      }
    }
    return true;
  }

private:
  const unsigned char* _pMat;
  size_t _nrow_fbm;
  std::vector<size_t> _rows;
  const int* _colInd;
  size_t _ncol;
  double _code256[256];
};

/******************************************************************************/

// A run of individuals (positions start to end-1 in GroupedRows::rows) that
// belong to the same group and use the same lookup table
struct RowRun {
  size_t start;
  size_t end;
  size_t group;
  size_t lut;
  bool contiguous; // TRUE if the rows are consecutive in the FBM
};

// Individuals reordered by group (and, within group, by lookup table), keeping
// their original relative order (stable sort), and split into runs
class GroupedRows {
public:
  // empty (to be assigned later)
  GroupedRows() {}

  // rows: 0-based FBM rows; group_ids and lut_ids: group and lookup table
  // of each individual; max_len: maximum length of a run (to avoid
  // overflowing the packed counters)
  GroupedRows(const std::vector<size_t>& rows0,
              const std::vector<size_t>& group_ids,
              const std::vector<size_t>& lut_ids,
              size_t max_len) {
    size_t n = rows0.size();
    // stable ordering by group and then lookup table
    order.resize(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
      if (group_ids[a] != group_ids[b]) return group_ids[a] < group_ids[b];
      return lut_ids[a] < lut_ids[b];
    });
    rows.resize(n);
    for (size_t k = 0; k < n; k++) rows[k] = rows0[order[k]];
    // split into runs of the same group and lookup table
    size_t k = 0;
    while (k < n) {
      size_t start = k;
      size_t g = group_ids[order[k]];
      size_t l = lut_ids[order[k]];
      while (k < n && group_ids[order[k]] == g && lut_ids[order[k]] == l &&
             (k - start) < max_len) {
        k++;
      }
      bool contiguous = true;
      for (size_t t = start + 1; t < k; t++) {
        if (rows[t] != rows[t - 1] + 1) contiguous = false;
      }
      runs.push_back({start, k, g, l, contiguous});
    }
  }

  std::vector<size_t> order; // order[k] = original position of k-th indiv
  std::vector<size_t> rows;  // 0-based FBM rows, sorted
  std::vector<RowRun> runs;
};

// convert 0-based group ids from R, checking that they are within range
inline std::vector<size_t> check_group_ids(const IntegerVector& groupIds,
                                           size_t n, size_t ngroups) {
  if (static_cast<size_t>(groupIds.size()) != n) {
    stop("groupIds must have one entry per individual (rowInd)");
  }
  std::vector<size_t> group_ids(n);
  for (size_t i = 0; i < n; i++) {
    if (groupIds[i] < 0 || static_cast<size_t>(groupIds[i]) >= ngroups) {
      stop("groupIds must be between 0 and ngroups - 1");
    }
    group_ids[i] = groupIds[i];
  }
  return group_ids;
}

/******************************************************************************/

// Sum packed lookup values over each run of individuals, for each locus.
//
// Each lookup table has 256 entries (one per byte value) of NFIELDS counters
// of 64 / NFIELDS bits each, packed in a 64-bit word (e.g. NFIELDS = 2 gives
// two 32-bit counters, NFIELDS = 3 three 21-bit counters). For locus j, the
// packed words of all individuals in a run are added up, and the counters
// are unpacked and added into tot[group * NFIELDS + field] (64-bit totals).
// Runs are short enough that the packed counters cannot overflow (see
// max_run_length(), to be used when building the GroupedRows).

// largest run of individuals for which the packed counters cannot overflow
template <int NFIELDS>
inline size_t max_run_length(const std::vector<uint64_t>& luts) {
  const int width = 64 / NFIELDS;
  const uint64_t mask = (width == 64) ? ~uint64_t(0)
                                      : ((uint64_t(1) << width) - 1);
  // largest value of any counter in any lookup table entry
  uint64_t max_value = 1;
  for (size_t e = 0; e < luts.size(); e++) {
    for (int k = 0; k < NFIELDS; k++) {
      max_value = std::max(max_value, (luts[e] >> (k * width)) & mask);
    }
  }
  return static_cast<size_t>(mask / max_value);
}

// Per-locus sums (pull interface): locus(j, tot) fills tot (ngroups * NFIELDS
// 64-bit totals) for locus j. It only reads plain C++ objects, so it can be
// called from multiple OpenMP threads, each with its own tot buffer. The
// accessor, grouped rows and lookup tables must outlive this object.
template <int NFIELDS>
class LutSums {
public:
  LutSums(const FBMByteAcc& acc,
          const GroupedRows& grouped,
          const std::vector<uint64_t>& luts,
          size_t ngroups)
    : _acc(acc), _grouped(grouped), _luts(luts), _ngroups(ngroups) {}

  size_t ngroups() const { return _ngroups; }

  inline void locus(size_t j, uint64_t* tot) const {
    const int width = 64 / NFIELDS;
    const uint64_t mask = (width == 64) ? ~uint64_t(0)
                                        : ((uint64_t(1) << width) - 1);
    const size_t* rows = _grouped.rows.data();
    const std::vector<RowRun>& runs = _grouped.runs;
    const unsigned char* col = _acc.col(j);
    std::fill(tot, tot + _ngroups * NFIELDS, 0);
    for (size_t r = 0; r < runs.size(); r++) {
      const RowRun& run = runs[r];
      const uint64_t* lut = _luts.data() + 256 * run.lut;
      uint64_t packed = 0;
      if (run.contiguous) {
        // consecutive rows: read the bytes directly
        const unsigned char* p = col + rows[run.start];
        size_t len = run.end - run.start;
        for (size_t i = 0; i < len; i++) packed += lut[p[i]];
      } else {
        for (size_t i = run.start; i < run.end; i++) {
          packed += lut[col[rows[i]]];
        }
      }
      // unpack the counters into the totals for this group
      uint64_t* t = tot + run.group * NFIELDS;
      for (int k = 0; k < NFIELDS; k++) {
        t[k] += (packed >> (k * width)) & mask;
      }
    }
  }

private:
  const FBMByteAcc& _acc;
  const GroupedRows& _grouped;
  const std::vector<uint64_t>& _luts;
  size_t _ngroups;
};

// Default block callback of lut_sums(), dosage_sums() and
// alt_valid_counts(), which does nothing.
struct NoBlockDone {
  void operator()(size_t, size_t) const {}
};

// Loop over all loci (push interface): for each locus j, compute the sums
// with LutSums::locus() and call finish(j, tot) to compute and store the
// statistics of that locus. Loci are processed in blocks of block_size,
// checking for user interrupts between blocks (outside of the parallel
// region). finish() is called in parallel, in no particular order; after
// each block, block_done(j0, j1) is called once (outside of the parallel
// region) with the range [j0, j1) of the loci of the block, so that results
// stored by finish() for that block can be combined in the order of the loci.
template <int NFIELDS, class Finish, class BlockDone = NoBlockDone>
void lut_sums(const FBMByteAcc& acc,
              const GroupedRows& grouped,
              const std::vector<uint64_t>& luts,
              size_t ngroups,
              int ncores,
              size_t block_size,
              Finish finish,
              BlockDone block_done = BlockDone()) {
  const LutSums<NFIELDS> sums(acc, grouped, luts, ngroups);
  const size_t m = acc.ncol();
  if (block_size < 1) block_size = 1;

  for (size_t j0 = 0; j0 < m; j0 += block_size) {
    size_t j1 = std::min(m, j0 + block_size);
#pragma omp parallel num_threads(ncores)
{
    // thread-local totals for one locus
    std::vector<uint64_t> tot(ngroups * NFIELDS);
#pragma omp for
    for (size_t j = j0; j < j1; j++) {
      sums.locus(j, tot.data());
      finish(j, tot.data());
    }
}
    block_done(j0, j1);
    Rcpp::checkUserInterrupt();
  }
}

/******************************************************************************/

// Weighted sums of decoded dosages (as doubles), for each locus and group.
//
// This is the fallback for codes that are not small integers (e.g. dosages),
// where integer lookup tables cannot be used. For each locus j and group g,
// tot[2 * g] is the sum of dosage * weight[i] and tot[2 * g + 1] the sum of
// copies[i] over non-missing genotypes (weight and copies are given in the
// original order of the individuals). Individuals are added in their original
// order within each group, so results are identical to accumulating over
// individuals in order. grouped should be built with all lut_ids equal.
// block_done(j0, j1) is called after each block, as in lut_sums().
template <class Finish, class BlockDone = NoBlockDone>
void dosage_sums(const FBMByteAcc& acc,
                 const GroupedRows& grouped,
                 const std::vector<double>& weight,
                 const std::vector<double>& copies,
                 size_t ngroups,
                 int ncores,
                 size_t block_size,
                 Finish finish,
                 BlockDone block_done = BlockDone()) {
  const size_t m = acc.ncol();
  const size_t n = grouped.rows.size();
  const size_t* rows = grouped.rows.data();
  const std::vector<RowRun>& runs = grouped.runs;
  // reorder the weights to match the sorted individuals
  std::vector<double> weight_sorted(n), copies_sorted(n);
  for (size_t k = 0; k < n; k++) {
    weight_sorted[k] = weight[grouped.order[k]];
    copies_sorted[k] = copies[grouped.order[k]];
  }
  double code256[256];
  for (int b = 0; b < 256; b++) code256[b] = acc.code(b);
  if (block_size < 1) block_size = 1;

  for (size_t j0 = 0; j0 < m; j0 += block_size) {
    size_t j1 = std::min(m, j0 + block_size);
#pragma omp parallel num_threads(ncores)
{
    std::vector<double> tot(ngroups * 2);
#pragma omp for
    for (size_t j = j0; j < j1; j++) {
      const unsigned char* col = acc.col(j);
      std::fill(tot.begin(), tot.end(), 0.0);
      for (size_t r = 0; r < runs.size(); r++) {
        const RowRun& run = runs[r];
        double* t = tot.data() + run.group * 2;
        for (size_t i = run.start; i < run.end; i++) {
          double x = code256[col[rows[i]]];
          if (x > -1) {
            t[0] += x * weight_sorted[i];
            t[1] += copies_sorted[i];
          }
        }
      }
      finish(j, tot.data());
    }
}
    block_done(j0, j1);
    Rcpp::checkUserInterrupt();
  }
}

/******************************************************************************/

// Weighted count of alternate alleles and number of valid allele copies, for
// each locus and group.
//
// For individual i with a non-missing genotype of dosage x, x * weight[i] is
// added to the alternate allele count and copies[i] (usually its ploidy) to
// the number of valid allele copies; weight[i] is 0.5 for pseudohaploids
// (whose single allele is stored as a 0/2 dosage) and 1 otherwise.
// finish(j, g, alt, valid) is called for each locus j and group g (in
// parallel), and block_done(j0, j1) after each block of loci (see
// lut_sums()).
//
// When all codes are integers and weights are 1 or 0.5, it uses one lookup
// table per (copies, weight) class, packing 2 * x * weight (an integer) and
// copies into two 32-bit counters: sums are then exact, so results are
// identical to summing the doubles. Otherwise (e.g. dosages from imputation
// probabilities) it falls back on summing doubles in the original order.
template <class Finish, class BlockDone = NoBlockDone>
void alt_valid_counts(const FBMByteAcc& acc,
                      const std::vector<size_t>& group_ids,
                      size_t ngroups,
                      const std::vector<double>& weight,
                      const std::vector<double>& copies,
                      int ncores,
                      size_t block_size,
                      Finish finish,
                      BlockDone block_done = BlockDone()) {
  const size_t n = acc.nrow();
  // check whether we can use integer lookup tables
  bool use_luts = acc.has_int_codes(255);
  for (size_t i = 0; i < n && use_luts; i++) {
    bool int_copies = copies[i] >= 0 && copies[i] <= 255 &&
      copies[i] == (int) copies[i];
    if (!int_copies || !(weight[i] == 1 || weight[i] == 0.5)) {
      use_luts = false;
    }
  }

  if (use_luts) {
    // find the distinct (copies, weight) classes, and the class of each
    // individual
    std::vector<double> class_copies, class_weight;
    std::vector<size_t> lut_ids(n);
    for (size_t i = 0; i < n; i++) {
      size_t c = 0;
      while (c < class_copies.size() &&
             !(class_copies[c] == copies[i] && class_weight[c] == weight[i])) {
        c++;
      }
      if (c == class_copies.size()) {
        class_copies.push_back(copies[i]);
        class_weight.push_back(weight[i]);
      }
      lut_ids[i] = c;
    }
    // one lookup table per class: low 32 bits hold 2 * x * weight, high 32
    // bits the number of copies; missing genotypes add nothing
    std::vector<uint64_t> luts(256 * class_copies.size(), 0);
    for (size_t c = 0; c < class_copies.size(); c++) {
      for (int b = 0; b < 256; b++) {
        double x = acc.code(b);
        if (x > -1) {
          uint64_t alt2 = static_cast<uint64_t>(2 * x * class_weight[c]);
          uint64_t n_copies = static_cast<uint64_t>(class_copies[c]);
          luts[256 * c + b] = alt2 | (n_copies << 32);
        }
      }
    }
    GroupedRows grouped(acc.rows(), group_ids, lut_ids,
                        max_run_length<2>(luts));
    lut_sums<2>(acc, grouped, luts, ngroups, ncores, block_size,
                [&](size_t j, const uint64_t* tot) {
                  for (size_t g = 0; g < ngroups; g++) {
                    finish(j, g, tot[2 * g] / 2.0,
                           static_cast<double>(tot[2 * g + 1]));
                  }
                }, block_done);
  } else {
    // sum doubles, keeping the original order of individuals within groups
    std::vector<size_t> lut_ids(n, 0);
    GroupedRows grouped(acc.rows(), group_ids, lut_ids,
                        std::numeric_limits<size_t>::max());
    dosage_sums(acc, grouped, weight, copies, ngroups, ncores, block_size,
                [&](size_t j, const double* tot) {
                  for (size_t g = 0; g < ngroups; g++) {
                    finish(j, g, tot[2 * g], tot[2 * g + 1]);
                  }
                }, block_done);
  }
}

/******************************************************************************/

// Summarise the genotypes of diploid (and pseudohaploid) individuals by group,
// one locus at a time (pull interface, for kernels that process loci in
// their own blocks, such as pop_stats_dip_pseudo_cpp() and
// pairwise_pop_fst_cpp()).
//
// For locus j (0-based index into colInd), locus() returns, for each group
// g, alt_sum[g], the sum of alternate allele dosages (multiplied by
// 1 / (3 - ploidy), i.e. halved for pseudohaploids), n_alleles[g], the
// number of valid alleles (sum of ploidy), and het2[g], twice the number of
// heterozygotes (dosage 1), over the individuals with a valid genotype.
//
// When all codes are integers and ploidy is 1 or 2, it uses one lookup
// table per ploidy that packs 2 * dosage * multiplier, ploidy and the
// heterozygote indicator into three 21-bit counters (LutSums<3>), so each
// genotype costs one integer addition; as the sums are exact, results are
// identical to summing doubles. Otherwise (e.g. dosages), it sums doubles by
// individual, in the original order within each group. locus() only reads
// plain C++ objects, so it can be called from multiple OpenMP threads, each
// with its own buffer of buffer_size() 64-bit integers.
class GroupedSummariesDipPseudo {
public:
  GroupedSummariesDipPseudo(Environment BM,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd,
                            const IntegerVector& groupIds,
                            size_t ngroups,
                            const NumericVector& ploidy)
    : _acc(BM, rowInd, colInd), _ngroups(ngroups) {
    size_t n = _acc.nrow();
    if (groupIds.size() != (R_xlen_t) n || ploidy.size() != (R_xlen_t) n) {
      stop("groupIds and ploidy should have one value per individual");
    }
    std::vector<size_t> group_ids(n);
    std::vector<double> multiplier(n), copies(n);
    for (size_t i = 0; i < n; i++) {
      if (groupIds[i] < 0 || (size_t) groupIds[i] >= ngroups) {
        stop("groupIds should be between 0 and ngroups - 1");
      }
      group_ids[i] = groupIds[i];
      copies[i] = ploidy[i];
      // multiplier to convert dosages for pseudohaploid data
      multiplier[i] = 1 / (3 - ploidy[i]);
    }

    // use lookup tables if codes are integers and ploidy is 1 or 2
    _use_luts = _acc.has_int_codes(255);
    for (size_t i = 0; i < n && _use_luts; i++) {
      if (!(copies[i] == 1 || copies[i] == 2)) _use_luts = false;
    }

    if (_use_luts) {
      // one lookup table per ploidy (table 0 for ploidy 1, 1 for ploidy 2):
      // 2 * dosage * multiplier, ploidy and the heterozygote indicator
      std::vector<size_t> lut_ids(n);
      for (size_t i = 0; i < n; i++) lut_ids[i] = (copies[i] == 1) ? 0 : 1;
      _luts.assign(256 * 2, 0);
      for (size_t c = 0; c < 2; c++) {
        double p = c + 1;
        for (int b = 0; b < 256; b++) {
          double x = _acc.code(b);
          // NA (and NaN) fail this comparison, so they are not valid
          if (x > -1) {
            uint64_t alt2 = static_cast<uint64_t>(2 * x * (1 / (3 - p)));
            uint64_t het = (x == 1) ? 1 : 0;
            _luts[256 * c + b] = alt2 | (static_cast<uint64_t>(p) << 21) |
              (het << 42);
          }
        }
      }
      _grouped = GroupedRows(_acc.rows(), group_ids, lut_ids,
                             max_run_length<3>(_luts));
    } else {
      // keep the original order of individuals within each group, with
      // their multiplier and ploidy, and decode bytes into doubles
      _grouped = GroupedRows(_acc.rows(), group_ids,
                             std::vector<size_t>(n, 0),
                             std::numeric_limits<size_t>::max());
      _multiplier.resize(n);
      _copies.resize(n);
      for (size_t k = 0; k < n; k++) {
        _multiplier[k] = multiplier[_grouped.order[k]];
        _copies[k] = copies[_grouped.order[k]];
      }
      for (int b = 0; b < 256; b++) {
        double x = _acc.code(b);
        _dosage[b] = x;
        _is_valid[b] = (x > -1);
        // we add 2 for heterozygotes so that we can then divide by valid
        // alleles (heterozygote counts are meaningless for pseudohaploids)
        _het_inc[b] = (x == 1) ? 2.0 : 0.0;
      }
    }
  }

  size_t nrow() const { return _acc.nrow(); }
  size_t ncol() const { return _acc.ncol(); }
  size_t ngroups() const { return _ngroups; }
  // size of the (per-thread) buffer needed by locus()
  size_t buffer_size() const { return 3 * _ngroups; }

  // Summarise locus j; each output array must have ngroups elements, and
  // buffer buffer_size() elements
  inline void locus(size_t j, double* alt_sum, double* n_alleles,
                    double* het2, uint64_t* buffer) const {
    if (_use_luts) {
      LutSums<3>(_acc, _grouped, _luts, _ngroups).locus(j, buffer);
      for (size_t g = 0; g < _ngroups; g++) {
        alt_sum[g] = buffer[3 * g] / 2.0;
        n_alleles[g] = static_cast<double>(buffer[3 * g + 1]);
        het2[g] = 2.0 * buffer[3 * g + 2];
      }
    } else {
      std::fill(alt_sum, alt_sum + _ngroups, 0.0);
      std::fill(n_alleles, n_alleles + _ngroups, 0.0);
      std::fill(het2, het2 + _ngroups, 0.0);
      const unsigned char* col = _acc.col(j);
      const std::vector<size_t>& rows = _grouped.rows;
      for (size_t r = 0; r < _grouped.runs.size(); r++) {
        const RowRun& run = _grouped.runs[r];
        size_t g = run.group;
        for (size_t k = run.start; k < run.end; k++) {
          unsigned char b = col[rows[k]];
          if (_is_valid[b]) {
            alt_sum[g] += _dosage[b] * _multiplier[k];
            n_alleles[g] += _copies[k];
            // adding 0 for homozygotes leaves the sum unchanged
            het2[g] += _het_inc[b];
          }
        }
      }
    }
  }

private:
  FBMByteAcc _acc;
  size_t _ngroups;
  bool _use_luts;
  GroupedRows _grouped;
  // lookup tables (integer codes)
  std::vector<uint64_t> _luts;
  // per-individual constants (in the order of _grouped.rows) and decoded
  // bytes (other codes)
  std::vector<double> _multiplier, _copies;
  double _dosage[256];
  bool _is_valid[256];
  double _het_inc[256];
};

/******************************************************************************/

#endif // TIDYPOPGEN_FBM_BYTE_ACC_H
