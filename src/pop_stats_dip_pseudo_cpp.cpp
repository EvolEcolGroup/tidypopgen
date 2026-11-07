/******************************************************************************/

#include "fbm_byte_acc.h"
#include <cmath>

/******************************************************************************/
// Population statistics for diploid (and pseudohaploid) data, following
// Nei 1987 (as in hierfstat::basic.stats()), computed in a single pass over
// the FBM.
//
// This is used by pop_global_stats(), pop_het_exp(), pop_het_obs() and
// pop_fis(). Loci are processed in blocks: for each block, we summarise the
// genotypes of each group, and compute (in parallel over loci) the
// statistics for each group (Ho, Hs, Fis) and the global statistics (Ho, Hs,
// Ht, Dst, Htp, Dstp, Fst, Fstp, Fis, Dest). Statistics by locus are only
// stored if requested; genome-wide means are accumulated locus by locus.
//
// The arithmetic follows, operation by operation, the R code that used to
// compute these statistics from loci x groups matrices, so results are
// identical: in particular, means across groups (rowMeans) and across loci
// (colMeans) are accumulated in long double, as R does.

// number of global statistics
const int n_global = 10;

// mean of the non-NaN values of x (length n, with stride), accumulated in
// long double as in R's rowMeans(x, na.rm = TRUE)
inline double mean_na_rm(const double* x, size_t n) {
  long double sum = 0.0;
  int cnt = 0;
  for (size_t i = 0; i < n; i++) {
    if (!std::isnan(x[i])) {
      sum += x[i];
      cnt++;
    }
  }
  sum /= cnt;
  return (double) sum;
}

// Running means across loci (as R's colMeans(x, na.rm = TRUE)). If
// skip_inf is true, infinite values are also ignored (as done by
// pop_global_stats(), which sets them to NA before taking the means)
class RunningMeans {
public:
  RunningMeans(size_t n, bool skip_inf) :
  _sum(n, 0.0), _cnt(n, 0), _skip_inf(skip_inf) {}

  inline void add(size_t k, double x) {
    if (!std::isnan(x) && !(_skip_inf && std::isinf(x))) {
      _sum[k] += x;
      _cnt[k]++;
    }
  }

  NumericVector means() const {
    NumericVector res(_sum.size());
    for (size_t k = 0; k < _sum.size(); k++) {
      long double sum = _sum[k];
      sum /= _cnt[k];
      res[k] = (double) sum;
    }
    return res;
  }

private:
  std::vector<long double> _sum;
  std::vector<R_xlen_t> _cnt;
  bool _skip_inf;
};

/******************************************************************************/

// [[Rcpp::export]]
List pop_stats_dip_pseudo_cpp(Environment BM,
                              const IntegerVector& rowInd,
                              const IntegerVector& colInd,
                              const IntegerVector& groupIds,
                              size_t ngroups,
                              const NumericVector& ploidy,
                              bool by_locus_ho,
                              bool by_locus_hs,
                              bool by_locus_fis,
                              bool by_locus_global,
                              int block_size,
                              int ncores) {

  // lookup tables and per-individual constants
  GroupedSummariesDipPseudo summ(BM, rowInd, colInd, groupIds, ngroups,
                                 ploidy);
  size_t m = summ.ncol(); // number of loci
  if (block_size < 1) stop("block_size should be a positive integer");

  // statistics by locus (empty matrices if not requested)
  NumericMatrix ho_locus(by_locus_ho ? m : 0, by_locus_ho ? ngroups : 0);
  NumericMatrix hs_locus(by_locus_hs ? m : 0, by_locus_hs ? ngroups : 0);
  NumericMatrix fis_locus(by_locus_fis ? m : 0, by_locus_fis ? ngroups : 0);
  NumericMatrix global_locus(by_locus_global ? m : 0,
                             by_locus_global ? n_global : 0);
  double* p_ho_locus = REAL(ho_locus);
  double* p_hs_locus = REAL(hs_locus);
  double* p_fis_locus = REAL(fis_locus);
  double* p_global_locus = REAL(global_locus);

  // genome-wide means of the group statistics (Ho, Hs, Fis for each group)
  RunningMeans group_means(3 * ngroups, false);
  // genome-wide means of the global statistics (keeping infinite values, as
  // when taking colMeans() of the statistics by locus, and dropping them, as
  // in pop_global_stats())
  RunningMeans global_means(n_global, false);
  RunningMeans global_means_finite(n_global, true);

  // statistics for the current block of loci; element [jj * ngroups + g] is
  // group g, locus jj of the block (and [jj * n_global + k] for the global
  // statistic k)
  size_t bs = std::min((size_t) block_size, std::max(m, (size_t) 1));
  std::vector<double> b_ho(bs * ngroups), b_hs(bs * ngroups),
    b_fis(bs * ngroups), b_global(bs * n_global);

  for (size_t j0 = 0; j0 < m; j0 += bs) {
    // number of loci in this block
    size_t nb = std::min(bs, m - j0);

#pragma omp parallel num_threads(ncores)
{
    // per-thread buffers for the group summaries
    std::vector<double> alt_sum(ngroups), n_alleles(ngroups), het2(ngroups);
    std::vector<uint64_t> summ_buffer(summ.buffer_size());
    std::vector<double> freq_alt(ngroups), freq_ref(ngroups), n_ind(ngroups);
    std::vector<double> sp2(ngroups);
    // buffers for the sums of the groups of four loci
    std::vector<double> alt_sum4(4 * ngroups), n_alleles4(4 * ngroups),
      het24(4 * ngroups);
    std::vector<uint64_t> summ_buffer4(summ.buffer_size4());
    // the statistics of locus j0 + jj, from the sums of its groups
    auto locus_stats = [&](size_t jj, const double* alt_sum,
                           const double* n_alleles, const double* het2) {
      size_t j = j0 + jj;
      double* ho = &b_ho[jj * ngroups];
      double* hs = &b_hs[jj * ngroups];
      double* fis = &b_fis[jj * ngroups];
      // statistics for each group
      for (size_t g = 0; g < ngroups; g++) {
        freq_alt[g] = alt_sum[g] / n_alleles[g];
        freq_ref[g] = 1 - freq_alt[g];
        ho[g] = het2[g] / n_alleles[g];
        // number of individuals
        n_ind[g] = n_alleles[g] / 2;
        // sum of squared frequencies
        sp2[g] = freq_alt[g] * freq_alt[g] + freq_ref[g] * freq_ref[g];
        // gene diversity (Nei 1987 eqn 7.39)
        hs[g] = (n_ind[g] / (n_ind[g] - 1)) *
          ((1 - sp2[g]) - ((ho[g] / 2) / n_ind[g]));
        fis[g] = 1 - ho[g] / hs[g];
      }
      // global statistics; means across groups are computed in long double
      // (as rowMeans() in R)
      double mHo = mean_na_rm(ho, ngroups);
      // number of groups with data (np) and harmonic mean of sample sizes
      // (mn); groups with no valid genotypes at this locus are treated as
      // missing (as hierfstat::ind.count returns NA for them)
      double np = 0.0, denom = 0.0;
      for (size_t g = 0; g < ngroups; g++) {
        if (!std::isnan(n_ind[g]) && n_ind[g] > 0.0) {
          np += 1.0;
          denom += 1.0 / n_ind[g];
        }
      }
      double mn = (denom > 0.0) ? np / denom : NA_REAL;
      // mean sum of squared frequencies, and sum of squared mean frequencies
      double msp2 = mean_na_rm(sp2.data(), ngroups);
      double mfreq_alt = mean_na_rm(freq_alt.data(), ngroups);
      double mfreq_ref = mean_na_rm(freq_ref.data(), ngroups);
      double mp2 = mfreq_alt * mfreq_alt + mfreq_ref * mfreq_ref;
      double mHs = (mn / (mn - 1)) * ((1 - msp2) - ((mHo / 2) / mn));
      double Ht = ((1 - mp2) + ((mHs / mn) / np)) - (((mHo / 2) / mn) / np);
      double mFis = 1 - mHo / mHs;
      double Dst = Ht - mHs;
      double Dstp = (np / (np - 1)) * Dst;
      double Htp = mHs + Dstp;
      double* glob = &b_global[jj * n_global];
      glob[0] = mHo;
      glob[1] = mHs;
      glob[2] = Ht;
      glob[3] = Dst;
      glob[4] = Htp;
      glob[5] = Dstp;
      glob[6] = Dst / Ht; // Fst
      glob[7] = Dstp / Htp; // Fstp
      glob[8] = mFis;
      glob[9] = Dstp / (1 - mHs); // Dest
      // store the statistics by locus, if requested
      for (size_t g = 0; g < ngroups; g++) {
        if (by_locus_ho) p_ho_locus[j + g * m] = ho[g];
        if (by_locus_hs) p_hs_locus[j + g * m] = hs[g];
        if (by_locus_fis) p_fis_locus[j + g * m] = fis[g];
      }
      if (by_locus_global) {
        for (int k = 0; k < n_global; k++) {
          p_global_locus[j + k * m] = glob[k];
        }
      }
    };
    // loci are summed four at a time, and the last (up to three) loci of the
    // block one at a time
    size_t n_quads = nb / 4;
#pragma omp for
    for (size_t q = 0; q < n_quads; q++) {
      summ.loci4(j0 + 4 * q, alt_sum4.data(), n_alleles4.data(),
                 het24.data(), summ_buffer4.data());
      for (size_t u = 0; u < 4; u++) {
        locus_stats(4 * q + u, &alt_sum4[u * ngroups],
                    &n_alleles4[u * ngroups], &het24[u * ngroups]);
      }
    }
#pragma omp for
    for (size_t jj = 4 * n_quads; jj < nb; jj++) {
      summ.locus(j0 + jj, alt_sum.data(), n_alleles.data(), het2.data(),
                 summ_buffer.data());
      locus_stats(jj, alt_sum.data(), n_alleles.data(), het2.data());
    }
}

    // update the genome-wide means, locus by locus (in order)
    for (size_t jj = 0; jj < nb; jj++) {
      for (size_t g = 0; g < ngroups; g++) {
        group_means.add(g, b_ho[jj * ngroups + g]);
        group_means.add(ngroups + g, b_hs[jj * ngroups + g]);
        group_means.add(2 * ngroups + g, b_fis[jj * ngroups + g]);
      }
      for (int k = 0; k < n_global; k++) {
        global_means.add(k, b_global[jj * n_global + k]);
        global_means_finite.add(k, b_global[jj * n_global + k]);
      }
    }
  }

  // split the means of the group statistics
  NumericVector all_group_means = group_means.means();
  NumericVector mean_ho(ngroups), mean_hs(ngroups), mean_fis(ngroups);
  for (size_t g = 0; g < ngroups; g++) {
    mean_ho[g] = all_group_means[g];
    mean_hs[g] = all_group_means[ngroups + g];
    mean_fis[g] = all_group_means[2 * ngroups + g];
  }

  return List::create(_["ho"] = ho_locus,
                      _["hs"] = hs_locus,
                      _["fis"] = fis_locus,
                      _["global"] = global_locus,
                      _["mean_ho"] = mean_ho,
                      _["mean_hs"] = mean_hs,
                      _["mean_fis"] = mean_fis,
                      _["mean_global"] = global_means.means(),
                      _["mean_global_finite"] = global_means_finite.means());
}

/******************************************************************************/
