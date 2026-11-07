/******************************************************************************/

#include "fbm_byte_acc.h"
#include <cmath>

/******************************************************************************/
// Runs of homozygosity (or heterozygosity) with sliding windows, as in
// detectRUNS::slidingRUNS.run() (and PLINK), used by windows_indiv_roh().
//
// For each individual, windows of window_size consecutive loci (in the order
// of the loci table) slide one locus at a time. A window is homozygous (for
// runs of homozygosity) if it has at most max_opp heterozygous genotypes, at
// most max_miss missing genotypes, and no gap between consecutive loci larger
// than max_gap. A locus is in a run if the proportion of homozygous windows
// among the windows that include it is larger than threshold. Runs are the
// stretches of consecutive loci in a run on each chromosome, which are then
// filtered by their number of loci, length, density and (optionally) their
// numbers of opposite and missing genotypes. For runs of heterozygosity, the
// roles of homozygous and heterozygous genotypes are swapped.
//
// detectRUNS works on one individual at a time, on a copy of its genotypes,
// and builds its runs with plyr; here, all individuals are processed in a
// single pass over the loci (reading the FBM column by column, as it is
// stored), keeping for each individual only the classes of its last
// window_size genotypes, the status of its last window_size windows, and its
// current run. The results are the same as detectRUNS, including features of
// its algorithm:
// - windows slide across chromosomes (only the gaps between consecutive loci
//   in the loci table are checked, and the gap between two chromosomes is
//   usually negative);
// - the proportion of windows is computed in single precision (float);
// - the number of windows that include each locus is that of detectRUNS,
//   which overestimates it for the first and last loci when there are fewer
//   than 2 * window_size - 1 loci;
// - genotypes are classified by their dosage: 0 and 2 are homozygous, 1 is
//   heterozygous, NA is missing, and other dosages are neither.
// With fewer loci than window_size, there are no windows, and no runs
// (detectRUNS fails). The loci of each chromosome must be adjacent in the
// loci table (checked in R), so that each individual has at most one open
// run at a time.

namespace {

// genotype classes
const unsigned char ROH_HOM = 0;
const unsigned char ROH_HET = 1;
const unsigned char ROH_MISS = 2;
const unsigned char ROH_OTHER = 3;

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
// and its numbers of loci, and of opposite and missing genotypes
struct OpenRun {
  bool open;
  uint32_t start;
  uint32_t last;
  int n_snp;
  int n_opp;
  int n_miss;
};

// a run that passed the filters (individual, chromosome code, first and last
// loci, number of loci)
struct FoundRun {
  uint32_t indiv;
  int32_t chrom;
  uint32_t start;
  uint32_t last;
  int32_t n_snp;
};

// parameters of the search, shared by all individuals
struct RohParams {
  size_t n_loci;
  size_t window_size;
  size_t n_windows;
  float threshold;
  unsigned char opp_class;
  int max_opp;
  int max_miss;
  double min_snp;
  double min_length_bps;
  double min_density;
  bool check_opp_run;
  double max_opp_run;
  bool check_miss_run;
  double max_miss_run;
};

// the number of windows that include locus k, as in detectRUNS
inline int roh_n_win(size_t k, size_t n_loci, size_t w) {
  bool first = (k < w);
  bool last = (n_loci - 1 - k < w);
  int a = (int) (k + 1);
  int b = (int) (n_loci - k);
  if (first && last) return std::max(a, b);
  if (first) return a;
  if (last) return b;
  return (int) w;
}

// what is known about the locus that is classified (finalised) at a step:
// its index, the slot of its genotype class (and of the status of the
// window that starts at it) in the rings of each individual, the smallest
// number of good windows that puts it in a run, its chromosome, and whether
// it is the last locus of its chromosome
struct LocusInfo {
  size_t k;
  size_t slot;
  int h_min;
  int chrom;
  bool chrom_end;
};

// the smallest number h of good windows (out of the n_win windows that
// include a locus) for which the locus is in a run, i.e. for which
// h / n_win > threshold, computed in single precision as in detectRUNS
// (h / n_win increases with h, so this is the same as testing each h); if
// no h is large enough, n_win + 1
inline int roh_h_min(int n_win, float threshold) {
  for (int h = 0; h <= n_win; h++) {
    float quotient = (float) h / n_win;
    if (quotient > threshold) return h;
  }
  return n_win + 1;
}

// the state of an individual (the classes of its last w loci and the status
// of its last w windows are in rings, in separate arrays)
struct RohState {
  int n_opp;   // opposite genotypes in the current window
  int n_miss;  // missing genotypes in the current window
  int h_sum;   // good windows among the last w windows
  OpenRun run; // the current run
};

// filter the run r of individual indiv, and keep it if it passes
inline void roh_close(OpenRun& r, size_t indiv, int chrom,
                      const RohParams& par, const Positions& pos,
                      std::vector<FoundRun>& found) {
  r.open = false;
  // as in detectRUNS::createRUNdf()
  double length_bps = pos[r.last] - pos[r.start];
  bool keep = r.n_snp >= par.min_snp &&
    length_bps >= par.min_length_bps &&
    (r.n_snp / length_bps) * 1000 >= par.min_density;
  if (par.check_opp_run && !(r.n_opp <= par.max_opp_run)) keep = false;
  if (par.check_miss_run && !(r.n_miss <= par.max_miss_run)) keep = false;
  if (keep) {
    FoundRun f = {(uint32_t) indiv, chrom, r.start, r.last, r.n_snp};
    found.push_back(f);
  }
}

// classify locus info.k (with genotype class cls), given h_sum, the number
// of good windows that include it, and update the run r
inline void roh_finalise(OpenRun& r, int h_sum, unsigned char cls,
                         const LocusInfo& info, size_t indiv,
                         const RohParams& par, const Positions& pos,
                         std::vector<FoundRun>& found) {
  if (h_sum >= info.h_min) {
    if (!r.open) {
      r.open = true;
      r.start = (uint32_t) info.k;
      r.n_snp = 0;
      r.n_opp = 0;
      r.n_miss = 0;
    }
    r.last = (uint32_t) info.k;
    r.n_snp++;
    r.n_opp += (cls == par.opp_class);
    r.n_miss += (cls == ROH_MISS);
    // the last locus of a chromosome ends its run
    if (info.chrom_end) roh_close(r, indiv, info.chrom, par, pos, found);
  } else if (r.open) {
    roh_close(r, indiv, info.chrom, par, pos, found);
  }
}

// Add loci t0 to t1 - 1 (of the block that starts at locus j0) for one
// individual (with its state st and its rings geno and status), reading its
// genotypes at offset row in the columns of the block (cols). For each locus
// t of the block, slots[t - j0] is its slot in the rings, and, once the
// window that ends at t is complete (gap_ok[t - j0] >= 0, whether that window
// has no large gap), infos[t - j0] describes the locus that is then
// finalised (the first locus of that window). The state is kept in local
// variables while the loci are added.
inline void roh_add_loci(RohState& st, unsigned char* geno,
                         unsigned char* status, size_t t0, size_t t1,
                         size_t j0, const unsigned char* const* cols,
                         size_t row, const unsigned char* lut,
                         const size_t* slots, const signed char* gap_ok,
                         const LocusInfo* infos, size_t indiv,
                         const RohParams& par, const Positions& pos,
                         std::vector<FoundRun>& found) {
  const size_t w = par.window_size;
  const unsigned char opp = par.opp_class;
  int n_opp = st.n_opp;
  int n_miss = st.n_miss;
  int h_sum = st.h_sum;
  OpenRun r = st.run;
  for (size_t t = t0; t < t1; t++) {
    size_t b = t - j0;
    unsigned char cls = lut[cols[b][row]];
    size_t slot_t = slots[b];
    // the locus that leaves the window (t - w) is replaced by locus t
    if (t >= w) {
      unsigned char old = geno[slot_t];
      n_opp -= (old == opp);
      n_miss -= (old == ROH_MISS);
    }
    geno[slot_t] = cls;
    n_opp += (cls == opp);
    n_miss += (cls == ROH_MISS);
    int window_gap_ok = gap_ok[b];
    if (window_gap_ok < 0) continue;
    // the window that starts at locus info.k is complete
    const LocusInfo& info = infos[b];
    unsigned char ok = (n_opp <= par.max_opp) & (n_miss <= par.max_miss) &
      (window_gap_ok == 1);
    if (info.k >= w) h_sum -= status[info.slot];
    status[info.slot] = ok;
    h_sum += ok;
    // locus info.k is in no later window, so it can now be classified
    roh_finalise(r, h_sum, geno[info.slot], info, indiv, par, pos, found);
  }
  st.n_opp = n_opp;
  st.n_miss = n_miss;
  st.h_sum = h_sum;
  st.run = r;
}

// after the last locus, finalise locus info.k for one individual: the window
// that starts at locus info.k - w (if any) no longer includes it
inline void roh_add_tail(RohState& st, const unsigned char* geno,
                         const unsigned char* status, const LocusInfo& info,
                         size_t indiv, const RohParams& par,
                         const Positions& pos, std::vector<FoundRun>& found) {
  if (info.k >= par.window_size) st.h_sum -= status[info.slot];
  roh_finalise(st.run, st.h_sum, geno[info.slot], info, indiv, par, pos,
               found);
}

} // namespace

/******************************************************************************/
// BM, rowInd, colInd: the genotypes (individuals as rows, loci in the order of
//   the loci table)
// chrom: the chromosome of each locus (1-based codes, e.g. of a factor); the
//   loci of each chromosome must be adjacent
// position: the position of each locus
// max_opp_run, max_miss_run: NA for no limit
// block_size: the number of loci read at a time (user interrupts are checked
//   between blocks)
// Returns a list with, for each run, the individual (1-based index into
// rowInd), the chromosome code, the first and last locus (1-based indices
// into colInd) and the number of loci in the run, ordered by individual and
// then along the loci table.

// [[Rcpp::export]]
List roh_sliding_cpp(Environment BM,
                     const IntegerVector& rowInd,
                     const IntegerVector& colInd,
                     const IntegerVector& chrom,
                     SEXP position,
                     int window_size,
                     double threshold,
                     bool heterozygosity,
                     int max_opp_window,
                     int max_miss_window,
                     double max_gap,
                     double min_snp,
                     double min_length_bps,
                     double min_density,
                     double max_opp_run,
                     double max_miss_run,
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
    } else if (g == 0 || g == 2) {
      lut[b] = ROH_HOM;
    } else if (g == 1) {
      lut[b] = ROH_HET;
    } else {
      lut[b] = ROH_OTHER;
    }
  }

  RohParams par;
  par.n_loci = m;
  par.window_size = window_size;
  par.n_windows = (m >= (size_t) window_size) ? m - window_size + 1 : 0;
  par.threshold = (float) threshold;
  par.opp_class = heterozygosity ? ROH_HOM : ROH_HET;
  par.max_opp = max_opp_window;
  par.max_miss = max_miss_window;
  par.min_snp = min_snp;
  par.min_length_bps = min_length_bps;
  par.min_density = min_density;
  par.check_opp_run = !ISNAN(max_opp_run);
  par.max_opp_run = max_opp_run;
  par.check_miss_run = !ISNAN(max_miss_run);
  par.max_miss_run = max_miss_run;
  const size_t w = window_size;

  // the last locus of each chromosome
  const int* p_chrom = chrom.begin();
  int n_chrom = 0;
  for (size_t j = 0; j < m; j++) {
    if (p_chrom[j] == NA_INTEGER || p_chrom[j] < 1) {
      stop("chrom should be positive integer codes");
    }
    n_chrom = std::max(n_chrom, p_chrom[j]);
  }
  std::vector<size_t> chrom_last(n_chrom, 0);
  for (size_t j = 0; j < m; j++) chrom_last[p_chrom[j] - 1] = j;

  // whether each window has no gap larger than max_gap between its loci (as
  // detectRUNS, gaps are truncated to integers, a missing gap is not large,
  // and windows span chromosomes); this is the same for all individuals
  Positions pos(position);
  double max_gap_trunc = std::trunc(max_gap);
  // is the gap between loci k and k + 1 large?
  auto large_gap = [&](size_t k) -> int {
    double gap = std::trunc(pos[k + 1] - pos[k]);
    return (!ISNAN(gap) && gap > max_gap_trunc) ? 1 : 0;
  };
  std::vector<signed char> window_gap_ok(par.n_windows, 1);
  if (par.n_windows > 0 && w >= 2) {
    // number of large gaps between the loci of window i (gaps i to
    // i + w - 2)
    int n_large = 0;
    for (size_t k = 0; k + 1 < w; k++) n_large += large_gap(k);
    for (size_t i = 0; i < par.n_windows; i++) {
      if (i > 0) {
        n_large += large_gap(i + w - 2) - large_gap(i - 1);
      }
      window_gap_ok[i] = (n_large == 0);
    }
  }

  // what is known about the locus finalised at each step (the same for all
  // individuals)
  auto locus_info = [&](size_t k) -> LocusInfo {
    LocusInfo info;
    info.k = k;
    info.slot = k % w;
    info.h_min = roh_h_min(roh_n_win(k, m, w), par.threshold);
    info.chrom = p_chrom[k];
    info.chrom_end = (chrom_last[p_chrom[k] - 1] == k);
    return info;
  };

  // individuals are processed in groups (in parallel), each with its own
  // list of runs; each individual has its state and its two rings
  const size_t group_size = 16;
  size_t n_groups = (n + group_size - 1) / group_size;
  std::vector<RohState> states(n);
  for (size_t i = 0; i < n; i++) {
    states[i].n_opp = 0;
    states[i].n_miss = 0;
    states[i].h_sum = 0;
    states[i].run.open = false;
  }
  std::vector<unsigned char> geno_rings(n * w, 0), status_rings(n * w, 0);
  std::vector<std::vector<FoundRun> > found(n_groups);
  const size_t* rows = acc.rows().data();
  // for each locus of a block: its column, its slot in the rings, whether
  // the window that ends at it is complete and has no large gap (-1 if it is
  // not complete), and the locus that is then finalised
  std::vector<const unsigned char*> cols;
  std::vector<size_t> slots;
  std::vector<signed char> gap_ok;
  std::vector<LocusInfo> infos;
  // loci are added in sub-blocks, so that the bytes of a sub-block for the
  // individuals of a group stay in cache while each individual (whose state
  // is kept in registers) goes through them
  const size_t sub_block = 256;

  for (size_t j0 = 0; j0 < m; j0 += block_size) {
    size_t j1 = std::min(m, j0 + (size_t) block_size);
    size_t nb = j1 - j0;
    cols.resize(nb);
    slots.resize(nb);
    gap_ok.resize(nb);
    infos.resize(nb);
    for (size_t t = j0; t < j1; t++) {
      cols[t - j0] = acc.col(t);
      slots[t - j0] = t % w;
      if (t + 1 >= w) {
        infos[t - j0] = locus_info(t + 1 - w);
        gap_ok[t - j0] = window_gap_ok[t + 1 - w];
      } else {
        gap_ok[t - j0] = -1;
      }
    }
#pragma omp parallel for num_threads(ncores) schedule(static)
    for (size_t g = 0; g < n_groups; g++) {
      size_t i0 = g * group_size;
      size_t i1 = std::min(n, i0 + group_size);
      for (size_t s0 = j0; s0 < j1; s0 += sub_block) {
        size_t s1 = std::min(j1, s0 + sub_block);
        for (size_t i = i0; i < i1; i++) {
          roh_add_loci(states[i], &geno_rings[i * w], &status_rings[i * w],
                       s0, s1, j0, cols.data(), rows[i], lut, slots.data(),
                       gap_ok.data(), infos.data(), i, par, pos, found[g]);
        }
      }
    }
    Rcpp::checkUserInterrupt();
  }
  // the last w - 1 loci of each individual
  if (m >= w) {
    infos.clear();
    for (size_t k = m + 1 - w; k < m; k++) infos.push_back(locus_info(k));
#pragma omp parallel for num_threads(ncores) schedule(static)
    for (size_t g = 0; g < n_groups; g++) {
      size_t i1 = std::min(n, (g + 1) * group_size);
      for (size_t i = g * group_size; i < i1; i++) {
        for (size_t k = 0; k < infos.size(); k++) {
          roh_add_tail(states[i], &geno_rings[i * w], &status_rings[i * w],
                       infos[k], i, par, pos, found[g]);
        }
      }
    }
  }

  // collect the runs, by individual
  size_t n_runs = 0;
  for (size_t g = 0; g < n_groups; g++) n_runs += found[g].size();
  IntegerVector r_indiv(n_runs), r_chrom(n_runs), r_start(n_runs),
    r_end(n_runs), r_n_snp(n_runs);
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
    }
    std::vector<FoundRun>().swap(found[g]);
  }
  return List::create(_["indiv"] = r_indiv, _["chrom"] = r_chrom,
                      _["start"] = r_start, _["end"] = r_end,
                      _["n_snp"] = r_n_snp);
}

/******************************************************************************/
