/******************************************************************************/

#include "fbm_byte_acc.h"
#include <cmath>

/******************************************************************************/
// Runs of homozygosity (or heterozygosity) with sliding windows, as in
// detectRUNS::slidingRUNS.run() (and PLINK), used by windows_indiv_roh().
//
// For each individual, windows of window_size consecutive loci slide one
// locus at a time along each chromosome (as in PLINK, windows restart at each
// chromosome, so no window holds the loci of two chromosomes). A window is homozygous (for
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
// current run. detectRUNS slides its windows along the whole loci table, so
// that near the end of a chromosome they also hold the first loci of the
// next one; here, the windows of each chromosome are those of detectRUNS run
// on that chromosome alone, so the results are those of detectRUNS run on
// each chromosome separately, including features of its algorithm:
// - the proportion of windows is computed in single precision (float);
// - the number of windows that include each locus is that of detectRUNS,
//   which overestimates it for the first and last loci of a chromosome with
//   fewer than 2 * window_size - 1 loci;
// - genotypes are classified by their dosage: 0 and 2 are homozygous, 1 is
//   heterozygous, NA is missing, and other dosages are neither.
// A chromosome with fewer loci than window_size has no windows, and no runs
// (detectRUNS fails). The loci of each chromosome must be adjacent in the
// loci table (checked in R), so that each individual has at most one open
// run at a time, and its state can be reset at the end of each chromosome.

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
// loci, number of loci, and of opposite and missing genotypes)
struct FoundRun {
  uint32_t indiv;
  int32_t chrom;
  uint32_t start;
  uint32_t last;
  int32_t n_snp;
  int32_t n_opp;
  int32_t n_miss;
};

// parameters of the search, shared by all individuals
struct RohParams {
  size_t window_size;
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

// the number of windows that include locus k (of n_loci on its chromosome),
// as in detectRUNS
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

// what is known about a locus that is classified (finalised): its index, the
// slot of its genotype class (and of the status of the window that starts at
// it) in the rings of each individual, whether a window (the one that starts
// window_size loci earlier on its chromosome) must be dropped from the count
// of good windows before it is classified, the smallest number of good
// windows that puts it in a run, its chromosome, and whether it is the last
// locus of its chromosome
struct LocusInfo {
  size_t k;
  size_t slot;
  bool drop;
  int h_min;
  int chrom;
  bool chrom_end;
};

// what is known about locus t when it is added: the slot of its genotype
// class in the rings, whether it replaces a locus (of the same chromosome)
// in the window, whether the window that ends at it is complete and has no
// large gap (-1 if it is not complete, i.e. it would start on the previous
// chromosome), and whether it is the last locus of its chromosome, followed
// by the loci of the chromosome that are then finalised (tails[tail_begin]
// to tails[tail_end - 1] of the block)
struct LocusStep {
  size_t slot;
  bool replaces;
  signed char gap_ok;
  bool chrom_end;
  size_t tail_begin;
  size_t tail_end;
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
    FoundRun f = {(uint32_t) indiv, chrom, r.start, r.last, r.n_snp,
                  r.n_opp, r.n_miss};
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
// genotypes at offset row in the columns of the block (cols). steps[t - j0]
// describes locus t, and, once the window that ends at t is complete,
// infos[t - j0] describes the locus that is then finalised (the first locus
// of that window). After the last locus of a chromosome, its last loci (in
// no complete window) are finalised, and the state is reset for the next
// chromosome. The state is kept in local variables while the loci are added.
inline void roh_add_loci(RohState& st, unsigned char* geno,
                         unsigned char* status, size_t t0, size_t t1,
                         size_t j0, const unsigned char* const* cols,
                         size_t row, const unsigned char* lut,
                         const LocusStep* steps, const LocusInfo* infos,
                         const LocusInfo* tails, size_t indiv,
                         const RohParams& par, const Positions& pos,
                         std::vector<FoundRun>& found) {
  const unsigned char opp = par.opp_class;
  int n_opp = st.n_opp;
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
      n_opp -= (old == opp);
      n_miss -= (old == ROH_MISS);
    }
    geno[step.slot] = cls;
    n_opp += (cls == opp);
    n_miss += (cls == ROH_MISS);
    if (step.gap_ok >= 0) {
      // the window that starts at locus info.k is complete
      const LocusInfo& info = infos[b];
      unsigned char ok = (n_opp <= par.max_opp) &
        (n_miss <= par.max_miss) & (step.gap_ok == 1);
      if (info.drop) h_sum -= status[info.slot];
      status[info.slot] = ok;
      h_sum += ok;
      // locus info.k is in no later window, so it can now be classified
      roh_finalise(r, h_sum, geno[info.slot], info, indiv, par, pos, found);
    }
    if (step.chrom_end) {
      // the last loci of the chromosome: for each, the window that starts
      // w loci earlier (if any) no longer includes it
      for (size_t k = step.tail_begin; k < step.tail_end; k++) {
        const LocusInfo& info = tails[k];
        if (info.drop) h_sum -= status[info.slot];
        roh_finalise(r, h_sum, geno[info.slot], info, indiv, par, pos,
                     found);
      }
      // the next chromosome starts with no window (its run, if any, was
      // closed by its last locus)
      n_opp = 0;
      n_miss = 0;
      h_sum = 0;
    }
  }
  st.n_opp = n_opp;
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
// position: the position of each locus
// max_opp_run, max_miss_run: NA for no limit
// block_size: the number of loci read at a time (user interrupts are checked
//   between blocks)
// Returns a list with, for each run, the individual (1-based index into
// rowInd), the chromosome code, the first and last locus (1-based indices
// into colInd), the number of loci in the run, and its numbers of opposite
// and missing genotypes, ordered by individual and then along the loci
// table.

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
  par.window_size = window_size;
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

  // the smallest number of good windows that puts a locus in a run, for each
  // number of windows that include it (at most 2 * w - 2, see roh_n_win())
  std::vector<int> h_min_table(2 * w + 1);
  for (size_t n_win = 0; n_win < h_min_table.size(); n_win++) {
    h_min_table[n_win] = roh_h_min((int) n_win, par.threshold);
  }

  // what is known about locus k when it is finalised (the same for all
  // individuals); ring slots and window counts are relative to the first
  // locus of its chromosome
  auto locus_info = [&](size_t k) -> LocusInfo {
    size_t c = p_chrom[k] - 1;
    size_t k_chrom = k - chrom_first[c];
    size_t n_chrom_loci = chrom_last[c] - chrom_first[c] + 1;
    LocusInfo info;
    info.k = k;
    info.slot = k_chrom % w;
    info.drop = (k_chrom >= w);
    info.h_min = h_min_table[roh_n_win(k_chrom, n_chrom_loci, w)];
    info.chrom = p_chrom[k];
    info.chrom_end = (chrom_last[c] == k);
    return info;
  };

  // whether each window has no gap larger than max_gap between its loci (as
  // in detectRUNS, gaps are truncated to integers); this is the same for all
  // individuals, and is computed as the loci are added, from the number of
  // large gaps among the last w - 1 gaps
  Positions pos(position);
  double max_gap_trunc = std::trunc(max_gap);
  // is the gap between loci k and k + 1 large?
  auto large_gap = [&](size_t k) -> int {
    double gap = std::trunc(pos[k + 1] - pos[k]);
    return (!ISNAN(gap) && gap > max_gap_trunc) ? 1 : 0;
  };
  int n_large = 0;

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
  // for each locus of a block: its column, what is known about it, and the
  // locus that is finalised once the window that ends at it is complete;
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
      // the large gaps of the window that ends at t (gaps t + 1 - w to
      // t - 1)
      if (t >= 1) n_large += large_gap(t - 1);
      if (t >= w) n_large -= large_gap(t - w);
      cols[b] = acc.col(t);
      LocusStep& step = steps[b];
      step.slot = t_chrom % w;
      step.replaces = (t_chrom >= w);
      if (t_chrom + 1 >= w) {
        infos[b] = locus_info(t + 1 - w);
        step.gap_ok = (n_large == 0);
      } else {
        step.gap_ok = -1;
      }
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
          roh_add_loci(states[i], &geno_rings[i * w], &status_rings[i * w],
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
    r_end(n_runs), r_n_snp(n_runs), r_n_opp(n_runs), r_n_miss(n_runs);
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
      r_n_opp[k] = f.n_opp;
      r_n_miss[k] = f.n_miss;
    }
    std::vector<FoundRun>().swap(found[g]);
  }
  return List::create(_["indiv"] = r_indiv, _["chrom"] = r_chrom,
                      _["start"] = r_start, _["end"] = r_end,
                      _["n_snp"] = r_n_snp, _["n_opp"] = r_n_opp,
                      _["n_miss"] = r_n_miss);
}

/******************************************************************************/
