/******************************************************************************/

#include "fbm_byte_acc.h"
#include <cmath>

/******************************************************************************/
// Runs of homozygosity with sliding windows, as PLINK 1.9 --homozyg
// (calc_homozyg() and roh_update() in plink_homozyg.c), used by
// windows_indiv_roh().
//
// For each individual, on each chromosome with N loci (N >= w), a window of
// w consecutive loci starts at each of the loci 0 to N - w (windows do not
// cross chromosomes, but do span gaps between loci). A window is a hit if it
// has at most max_het_window heterozygous and at most max_miss_window missing
// genotypes. Each locus c is then classified as in a run if at least h_min of
// the windows that include it are hits, where h_min is
// (int) (n * threshold + 1 - EPSILON) (i.e. the ceiling of n * threshold), and
// n is min(c + 1, w) for c <= N - w, and N - c for the last w - 1 loci, as in
// PLINK. The loci are then scanned in order, and a run (if open) is closed
// before locus c if c is not in a run, if c is heterozygous and the run
// already has max_het heterozygous genotypes, or if the gap between loci c - 1
// and c is larger than max_gap; a run is opened at each locus in a run when
// none is open (except at a heterozygous locus if max_het is 0). A closed run
// is kept if it has at least min_snp loci, a length (last position - first
// position + 1, or without the + 1 with old lengths) of at least min_bases,
// and at most max_bases_per_snp bases per locus. The last run of each
// chromosome is closed at its last locus. A chromosome with fewer than w loci
// has no runs. Genotypes are classified by their dosage: 1 is heterozygous,
// NA is missing, and anything else is homozygous.
//
// All individuals are processed in a single pass over the loci (reading the
// FBM column by column, as it is stored), keeping for each individual only
// the classes of its last w genotypes, the hits of its last w windows, and
// its current run. The loci of each chromosome must be adjacent in the loci
// table (checked in R), so that the state of each individual can be reset at
// the end of each chromosome.

namespace {

// PLINK's tolerance in floating point comparisons (2^-30)
const double ROH_EPSILON = 0.000000000931322574615478515625;

// genotype classes
const unsigned char ROH_HOM = 0;
const unsigned char ROH_HET = 1;
const unsigned char ROH_MISS = 2;

// the positions of the loci, read from an integer or a double vector without
// copying it
struct Positions {
  const int* p_int;
  const double* p_double;
  Positions(SEXP x) : p_int(NULL), p_double(NULL) {
    if (TYPEOF(x) == INTSXP) {
      p_int = INTEGER(x);
    } else if (TYPEOF(x) == REALSXP) {
      p_double = REAL(x);
    } else {
      stop("position should be an integer or a numeric vector");
    }
  }
  inline double operator[](size_t k) const {
    if (p_int != NULL) {
      return (p_int[k] == NA_INTEGER) ? NA_REAL : (double) p_int[k];
    }
    return p_double[k];
  }
};

// the current run of an individual (if open), with its first and last loci,
// and its numbers of loci, and of heterozygous and missing genotypes
struct OpenRun {
  bool open;
  uint32_t start;
  uint32_t last;
  int n_snp;
  int n_het;
  int n_miss;
};

// a run that passed the filters (individual, chromosome code, first and last
// loci, number of loci, and of heterozygous and missing genotypes)
struct FoundRun {
  uint32_t indiv;
  int32_t chrom;
  uint32_t start;
  uint32_t last;
  int32_t n_snp;
  int32_t n_het;
  int32_t n_miss;
};

// parameters of the search, shared by all individuals
struct RohParams {
  size_t window_size;
  int max_het_window;
  int max_miss_window;
  double min_snp;
  double min_bases;
  double max_bases_per_snp;
  double max_het; // the largest number of heterozygous genotypes in a run
  double length_offset; // 1, or 0 for old lengths
};

// what is known about a locus that is classified: its index, the slot of
// its genotype class (and of the hit of the window that starts at it) in the
// rings of each individual, whether a window (the one that starts
// window_size loci earlier on its chromosome) must be dropped from the count
// of hits before it is classified, the smallest number of hits that puts it
// in a run, whether the gap from the previous locus is larger than max_gap,
// its chromosome, and whether it is the last locus of its chromosome
struct LocusInfo {
  size_t k;
  size_t slot;
  bool drop;
  int h_min;
  bool gap_before;
  int chrom;
  bool chrom_end;
};

// what is known about locus t when it is added: the slot of its genotype
// class in the rings, whether it replaces a locus (of the same chromosome)
// in the window, whether the window that ends at it is complete (so that the
// locus that starts it can be classified), and whether it is the last locus
// of its chromosome, followed by the loci of the chromosome that are then
// classified (tails[tail_begin] to tails[tail_end - 1] of the block)
struct LocusStep {
  size_t slot;
  bool replaces;
  bool complete;
  bool chrom_end;
  size_t tail_begin;
  size_t tail_end;
};

// the state of an individual (the classes of its last w loci and the hits
// of its last w windows are in rings, in separate arrays)
struct RohState {
  int n_het;   // heterozygous genotypes in the current window
  int n_miss;  // missing genotypes in the current window
  int h_sum;   // hits among the windows that include the next locus
  OpenRun run; // the current run
};

// close the run r of individual indiv, and keep it if it passes the filters
inline void roh_close(OpenRun& r, size_t indiv, int chrom,
                      const RohParams& par, const Positions& pos,
                      std::vector<FoundRun>& found) {
  r.open = false;
  double base_len = pos[r.last] + par.length_offset - pos[r.start];
  bool keep = r.n_snp >= par.min_snp && base_len >= par.min_bases &&
    (double) r.n_snp * par.max_bases_per_snp >= base_len;
  if (keep) {
    FoundRun f = {(uint32_t) indiv, chrom, r.start, r.last, r.n_snp,
                  r.n_het, r.n_miss};
    found.push_back(f);
  }
}

// classify locus info.k (with genotype class cls), given h_sum, the number
// of hits among the windows that include it, and update the run r, as
// roh_update() in PLINK
inline void roh_finalise(OpenRun& r, int h_sum, unsigned char cls,
                         const LocusInfo& info, size_t indiv,
                         const RohParams& par, const Positions& pos,
                         std::vector<FoundRun>& found) {
  bool hit = h_sum >= info.h_min;
  bool het = (cls == ROH_HET);
  if (r.open && (!hit || (het && r.n_het == par.max_het) ||
                 info.gap_before)) {
    roh_close(r, indiv, info.chrom, par, pos, found);
  }
  if (hit) {
    if (!r.open) {
      // a run can not start with a heterozygous genotype if none is allowed
      if (het && par.max_het == 0) return;
      r.open = true;
      r.start = (uint32_t) info.k;
      r.n_snp = 0;
      r.n_het = 0;
      r.n_miss = 0;
    }
    r.last = (uint32_t) info.k;
    r.n_snp++;
    r.n_het += het;
    r.n_miss += (cls == ROH_MISS);
  }
  // the last locus of a chromosome ends its run
  if (info.chrom_end && r.open) roh_close(r, indiv, info.chrom, par, pos, found);
}

// Add loci t0 to t1 - 1 (of the block that starts at locus j0) for one
// individual (with its state st and its rings geno and hits), reading its
// genotypes at offset row in the columns of the block (cols). steps[t - j0]
// describes locus t, and, once the window that ends at t is complete,
// infos[t - j0] describes the locus that is then classified (the first locus
// of that window). After the last locus of a chromosome, its last loci (in
// no complete window) are classified, and the state is reset for the next
// chromosome. The state is kept in local variables while the loci are added.
inline void roh_add_loci(RohState& st, unsigned char* geno,
                         unsigned char* hits, size_t t0, size_t t1,
                         size_t j0, const unsigned char* const* cols,
                         size_t row, const unsigned char* lut,
                         const LocusStep* steps, const LocusInfo* infos,
                         const LocusInfo* tails, size_t indiv,
                         const RohParams& par, const Positions& pos,
                         std::vector<FoundRun>& found) {
  int n_het = st.n_het;
  int n_miss = st.n_miss;
  int h_sum = st.h_sum;
  OpenRun r = st.run;
  for (size_t t = t0; t < t1; t++) {
    size_t b = t - j0;
    const LocusStep& step = steps[b];
    unsigned char cls = lut[cols[b][row]];
    // the locus that leaves the window (w loci earlier on the chromosome) is
    // replaced by locus t
    if (step.replaces) {
      unsigned char old = geno[step.slot];
      n_het -= (old == ROH_HET);
      n_miss -= (old == ROH_MISS);
    }
    geno[step.slot] = cls;
    n_het += (cls == ROH_HET);
    n_miss += (cls == ROH_MISS);
    if (step.complete) {
      // the window that starts at locus info.k is complete
      const LocusInfo& info = infos[b];
      unsigned char hit = (n_het <= par.max_het_window) &
        (n_miss <= par.max_miss_window);
      if (info.drop) h_sum -= hits[info.slot];
      hits[info.slot] = hit;
      h_sum += hit;
      // locus info.k is in no later window, so it can now be classified
      roh_finalise(r, h_sum, geno[info.slot], info, indiv, par, pos, found);
    }
    if (step.chrom_end) {
      // the last loci of the chromosome: for each, the window that starts
      // w loci earlier (if any) no longer includes it
      for (size_t k = step.tail_begin; k < step.tail_end; k++) {
        const LocusInfo& info = tails[k];
        if (info.drop) h_sum -= hits[info.slot];
        roh_finalise(r, h_sum, geno[info.slot], info, indiv, par, pos,
                     found);
      }
      // the next chromosome starts with no window (its run, if any, was
      // closed by its last locus)
      n_het = 0;
      n_miss = 0;
      h_sum = 0;
    }
  }
  st.n_het = n_het;
  st.n_miss = n_miss;
  st.h_sum = h_sum;
  st.run = r;
}

} // namespace

/******************************************************************************/
// BM, rowInd, colInd: the genotypes (individuals as rows, loci in the order of
//   the loci table)
// chrom: the chromosome of each locus (1-based codes, e.g. of a factor); the
//   loci of each chromosome must be adjacent
// position: the position of each locus (in bp)
// window_size, max_het_window, max_miss_window, threshold: the windows
// max_gap: the largest gap (in bp) between consecutive loci of a run
// min_snp, min_bases, max_bases_per_snp: filters of the runs
// max_het: the largest number of heterozygous genotypes in a run (Inf for no
//   limit)
// old_lengths: whether the length of a run is last - first position (rather
//   than last - first + 1)
// block_size: the number of loci read at a time (user interrupts are checked
//   between blocks)
// Returns a list with, for each run, the individual (1-based index into
// rowInd), the chromosome code, the first and last locus (1-based indices
// into colInd), the number of loci in the run, and its numbers of
// heterozygous and missing genotypes, ordered by individual and then along
// the loci table.

// [[Rcpp::export]]
List roh_sliding_cpp(Environment BM,
                     const IntegerVector& rowInd,
                     const IntegerVector& colInd,
                     const IntegerVector& chrom,
                     SEXP position,
                     int window_size,
                     int max_het_window,
                     int max_miss_window,
                     double threshold,
                     double max_gap,
                     double min_snp,
                     double min_bases,
                     double max_bases_per_snp,
                     double max_het,
                     bool old_lengths,
                     int ncores,
                     int block_size) {
  FBMByteAcc acc(BM, rowInd, colInd);
  size_t n = acc.nrow();
  size_t m = acc.ncol();
  if ((size_t) chrom.size() != m || (size_t) Rf_xlength(position) != m) {
    stop("chrom and position should have one value per locus");
  }
  if (m >= (size_t) std::numeric_limits<uint32_t>::max() ||
      n >= (size_t) std::numeric_limits<uint32_t>::max()) {
    stop("too many loci or individuals");
  }
  if (window_size < 1) stop("window_size should be at least 1");
  if (block_size < 1) block_size = 1;

  // the class of each byte of the FBM (through its code256)
  unsigned char lut[256];
  for (int b = 0; b < 256; b++) {
    double g = acc.code(b);
    if (ISNAN(g)) {
      lut[b] = ROH_MISS;
    } else if (g == 1) {
      lut[b] = ROH_HET;
    } else {
      lut[b] = ROH_HOM;
    }
  }

  RohParams par;
  par.window_size = window_size;
  par.max_het_window = max_het_window;
  par.max_miss_window = max_miss_window;
  par.min_snp = min_snp;
  par.min_bases = min_bases;
  par.max_bases_per_snp = max_bases_per_snp;
  par.max_het = max_het;
  par.length_offset = old_lengths ? 0 : 1;
  const size_t w = window_size;

  // the first and last locus of each chromosome
  const int* p_chrom = chrom.begin();
  int n_chrom = 0;
  for (size_t j = 0; j < m; j++) {
    if (p_chrom[j] == NA_INTEGER || p_chrom[j] < 1) {
      stop("chrom should be positive integer codes");
    }
    n_chrom = std::max(n_chrom, p_chrom[j]);
  }
  std::vector<size_t> chrom_first(n_chrom, 0), chrom_last(n_chrom, 0);
  for (size_t j = m; j-- > 0;) chrom_first[p_chrom[j] - 1] = j;
  for (size_t j = 0; j < m; j++) chrom_last[p_chrom[j] - 1] = j;

  // what is known about locus k when it is classified (the same for all
  // individuals); ring slots and window counts are relative to the first
  // locus of its chromosome
  Positions pos(position);
  auto locus_info = [&](size_t k) -> LocusInfo {
    size_t c = p_chrom[k] - 1;
    size_t k_chrom = k - chrom_first[c];
    size_t n_chrom_loci = chrom_last[c] - chrom_first[c] + 1;
    LocusInfo info;
    info.k = k;
    info.slot = k_chrom % w;
    info.drop = (k_chrom >= w);
    // the number of windows used for the threshold, as in PLINK (for the
    // last w - 1 loci, it is that of the loci to the end of the chromosome,
    // which is more than the windows that include the locus on chromosomes
    // with fewer than 2 * w - 2 loci)
    size_t n_win = (k_chrom + w <= n_chrom_loci) ?
      std::min(k_chrom + 1, w) : n_chrom_loci - k_chrom;
    info.h_min = (int) ((double) n_win * threshold + 1.0 - ROH_EPSILON);
    info.gap_before = k_chrom > 0 && (pos[k] - pos[k - 1]) > max_gap;
    info.chrom = p_chrom[k];
    info.chrom_end = (chrom_last[c] == k);
    return info;
  };

  // individuals are processed in groups (in parallel), each with its own
  // list of runs; each individual has its state and its two rings
  const size_t group_size = 16;
  size_t n_groups = (n + group_size - 1) / group_size;
  std::vector<RohState> states(n);
  for (size_t i = 0; i < n; i++) {
    states[i].n_het = 0;
    states[i].n_miss = 0;
    states[i].h_sum = 0;
    states[i].run.open = false;
  }
  std::vector<unsigned char> geno_rings(n * w, 0), hit_rings(n * w, 0);
  std::vector<std::vector<FoundRun> > found(n_groups);
  const size_t* rows = acc.rows().data();
  // for each locus of a block: its column, what is known about it, and the
  // locus that is classified once the window that ends at it is complete;
  // and the last loci of the chromosomes that end in the block
  std::vector<const unsigned char*> cols;
  std::vector<LocusStep> steps;
  std::vector<LocusInfo> infos;
  std::vector<LocusInfo> tails;
  // loci are added in sub-blocks, so that the bytes of a sub-block for the
  // individuals of a group stay in cache while each individual (whose state
  // is kept in registers) goes through them
  const size_t sub_block = 256;

  for (size_t j0 = 0; j0 < m; j0 += block_size) {
    size_t j1 = std::min(m, j0 + (size_t) block_size);
    size_t nb = j1 - j0;
    cols.resize(nb);
    steps.resize(nb);
    infos.resize(nb);
    tails.clear();
    for (size_t t = j0; t < j1; t++) {
      size_t b = t - j0;
      size_t c = p_chrom[t] - 1;
      size_t t_chrom = t - chrom_first[c];
      cols[b] = acc.col(t);
      LocusStep& step = steps[b];
      step.slot = t_chrom % w;
      step.replaces = (t_chrom >= w);
      step.complete = (t_chrom + 1 >= w);
      if (step.complete) infos[b] = locus_info(t + 1 - w);
      step.chrom_end = (chrom_last[c] == t);
      step.tail_begin = tails.size();
      // the last w - 1 loci of the chromosome, if it has any window
      if (step.chrom_end && t_chrom + 1 >= w) {
        for (size_t k = t + 2 - w; k <= t; k++) {
          tails.push_back(locus_info(k));
        }
      }
      step.tail_end = tails.size();
    }
#pragma omp parallel for num_threads(ncores) schedule(static)
    for (size_t g = 0; g < n_groups; g++) {
      size_t i0 = g * group_size;
      size_t i1 = std::min(n, i0 + group_size);
      for (size_t s0 = j0; s0 < j1; s0 += sub_block) {
        size_t s1 = std::min(j1, s0 + sub_block);
        for (size_t i = i0; i < i1; i++) {
          roh_add_loci(states[i], &geno_rings[i * w], &hit_rings[i * w],
                       s0, s1, j0, cols.data(), rows[i], lut, steps.data(),
                       infos.data(), tails.data(), i, par, pos, found[g]);
        }
      }
    }
    Rcpp::checkUserInterrupt();
  }
  // collect the runs, by individual
  size_t n_runs = 0;
  for (size_t g = 0; g < n_groups; g++) n_runs += found[g].size();
  IntegerVector r_indiv(n_runs), r_chrom(n_runs), r_start(n_runs),
    r_end(n_runs), r_n_snp(n_runs), r_n_het(n_runs), r_n_miss(n_runs);
  size_t k = 0;
  for (size_t g = 0; g < n_groups; g++) {
    // runs of the same group are ordered by end, not by individual
    std::stable_sort(found[g].begin(), found[g].end(),
                     [](const FoundRun& a, const FoundRun& b) {
                       return a.indiv < b.indiv;
                     });
    for (size_t r = 0; r < found[g].size(); r++, k++) {
      const FoundRun& f = found[g][r];
      r_indiv[k] = (int) f.indiv + 1;
      r_chrom[k] = f.chrom;
      r_start[k] = (int) f.start + 1;
      r_end[k] = (int) f.last + 1;
      r_n_snp[k] = f.n_snp;
      r_n_het[k] = f.n_het;
      r_n_miss[k] = f.n_miss;
    }
    std::vector<FoundRun>().swap(found[g]);
  }
  return List::create(_["indiv"] = r_indiv, _["chrom"] = r_chrom,
                      _["start"] = r_start, _["end"] = r_end,
                      _["n_snp"] = r_n_snp, _["n_het"] = r_n_het,
                      _["n_miss"] = r_n_miss);
}

/******************************************************************************/
