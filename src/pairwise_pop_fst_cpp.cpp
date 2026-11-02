/******************************************************************************/

#include "grouped_summaries_dip_pseudo.h"
#include <cmath>
#include <string>

/******************************************************************************/
// Numerator and denominator of pairwise Fst at a single locus, for each
// method. The arithmetic follows, operation by operation, the previous
// vectorised implementations (one loop over all loci per pair of
// populations), so that results are unchanged.

// Hudson's formulation (Bhatia et al 2013); h1 and h2 are the precomputed
// within population terms (freq_alt * freq_ref) / (n - 1)
inline void fst_hudson_num_den(double fa1, double fr1, double h1,
                               double fa2, double fr2, double h2,
                               double& num, double& den) {
  double diff = fa1 - fa2;
  num = diff * diff - h1 - h2;
  den = fa1 * fr2 + fa2 * fr1;
}

// Nei 1987 (as in hierfstat::pairwise.neifst()); n is the number of valid
// alleles, ho the observed heterozygosity
inline void fst_nei87_num_den(double n1, double ho1, double fa1, double fr1,
                              double n2, double ho2, double fa2, double fr2,
                              double& num, double& den) {
  // number of individuals
  double n_pair[2] = {n1 / 2.0, n2 / 2.0};
  double sho[2] = {ho1, ho2};
  // sample size and het_obs means
  int valid = 0;
  double nsum = 0.0, inv_nsum = 0.0;
  double ho_sum = 0.0;
  for (int j = 0; j < 2; ++j) {
    if (!std::isnan(n_pair[j])) {
      valid++;
      ho_sum += sho[j];
      nsum += 1.0;
      inv_nsum += 1.0 / n_pair[j];
    }
  }
  double np = valid;
  double mn = (inv_nsum > 0) ? nsum / inv_nsum : NA_REAL;
  double mHo = ho_sum / valid;
  // mean sum of squared frequencies
  double sp2a = fa1 * fa1 + fa2 * fa2;
  double sp2r = fr1 * fr1 + fr2 * fr2;
  double sp2 = sp2a + sp2r;
  double msp2 = sp2 / 2.0;
  // squared mean frequencies
  double freqAm = (fa1 + fa2) / 2.0;
  double freqRm = (fr1 + fr2) / 2.0;
  double mp2 = freqAm * freqAm + freqRm * freqRm;
  double mHs = mn / (mn - 1.0) * (1.0 - msp2 - mHo / (2.0 * mn));
  double Ht = 1.0 - mp2 + mHs / (mn * np) - mHo / (2.0 * mn * np);
  double Dst = Ht - mHs;
  double Dstp = (np / (np - 1.0)) * Dst;
  double Htp = mHs + Dstp;
  num = Dstp;
  den = Htp;
}

// Weir and Cockerham 1984 (as in hierfstat::pairwise.WCfst()); an is the
// number of valid alleles, p the alternative allele frequency, and h the
// observed heterozygosity
inline void fst_wc84_num_den(double an1, double p1, double h1,
                             double an2, double p2, double h2,
                             double& num, double& den) {
  const int r = 2; // number of populations
  double n_ind[2] = {an1 / 2.0, an2 / 2.0};
  double p[2] = {p1, p2};
  double h[2] = {h1, h2};
  double sum_n = 0.0, sum_sq = 0.0;
  for (int j = 0; j < r; ++j) {
    sum_n += n_ind[j];
    sum_sq += n_ind[j] * n_ind[j];
  }
  double n_total = sum_n;
  double n_bar = sum_n / r;
  double n_c = (sum_n - sum_sq / sum_n) / (r - 1);
  double sum_pn = 0.0, sum_sq_diff = 0.0, sum_h = 0.0;
  for (int j = 0; j < r; ++j) {
    sum_pn += p[j] * n_ind[j];
    sum_h += h[j] * n_ind[j];
  }
  double p_bar = sum_pn / n_total;
  double h_bar = sum_h / n_total;
  for (int j = 0; j < r; ++j) {
    double diff = p[j] - p_bar;
    sum_sq_diff += diff * diff * n_ind[j];
  }
  double s2 = sum_sq_diff / (n_bar * (r - 1));
  double a = n_bar / n_c * (s2 -
    (1.0 / (n_bar - 1.0)) *
    (p_bar * (1 - p_bar) - ((r - 1.0) / r) * s2 - h_bar / 4.0));
  double b = n_bar / (n_bar - 1.0) *
    (p_bar * (1 - p_bar) - ((r - 1.0) / r) * s2 -
    ((2 * n_bar - 1.0) / (4.0 * n_bar)) * h_bar);
  double c = h_bar / 2.0;
  num = a;
  den = a + b + c;
}

/******************************************************************************/
// Helpers to check the arguments shared by the exported functions

// convert the name of the method to an integer id
inline int fst_method_id(std::string method) {
  if (method == "Hudson") {
    return 0;
  } else if (method == "Nei87") {
    return 1;
  } else if (method == "WC84") {
    return 2;
  }
  stop("method should be one of 'Hudson', 'Nei87' or 'WC84'");
  return -1; // never reached
}

// convert the pairs of groups to 0-based indices, checking them
inline void fst_pairs(const IntegerMatrix& pairwise_combn, size_t ngroups,
                      std::vector<int>& pop1, std::vector<int>& pop2) {
  size_t n_combn = pairwise_combn.ncol();
  pop1.resize(n_combn);
  pop2.resize(n_combn);
  for (size_t k = 0; k < n_combn; k++) {
    pop1[k] = pairwise_combn(0, k) - 1;
    pop2[k] = pairwise_combn(1, k) - 1;
    if (pop1[k] < 0 || pop2[k] < 0 || (size_t) pop1[k] >= ngroups ||
        (size_t) pop2[k] >= ngroups) {
      stop("pairwise_combn should contain group indices between 1 and ngroups");
    }
  }
}

/******************************************************************************/
// Compute the numerator and denominator of pairwise Fst for all pairs of
// groups, in a single pass over the FBM.
//
// Loci are processed in blocks of `block_size` loci: for each block, we first
// summarise the genotypes of each group (in parallel over loci), and then
// compute the numerator and denominator of Fst for every pair of groups (in
// parallel over pairs). For each pair, consume(k, j0, nb, num, den) is then
// called with the numerators and denominators of the nb loci of the block
// (starting at locus j0); as each pair is handled by a single thread, and
// blocks are processed in order, consume() sees the loci of each pair in
// order, so sums accumulated by it do not depend on the number of cores or
// on the block size. Memory use is limited to the block summaries, as we
// never store the full loci x groups summaries.
template <class Consumer>
void pairwise_fst_blocks(const GroupedSummariesDipPseudo& summ,
                         const std::vector<int>& pop1,
                         const std::vector<int>& pop2,
                         int method_id,
                         int block_size,
                         int ncores,
                         Consumer consume) {
  if (block_size < 1) stop("block_size should be a positive integer");
  size_t m = summ.ncol(); // number of loci
  size_t ngroups = summ.ngroups();
  size_t n_combn = pop1.size();

  // summaries of the current block of loci; each array is stored by group,
  // with the loci of the block contiguous for each group (i.e. element
  // [g * bs + jj] is group g, locus jj of the block)
  size_t bs = std::min((size_t) block_size, std::max(m, (size_t) 1));
  std::vector<double> b_freq_alt(ngroups * bs), b_freq_ref(ngroups * bs);
  std::vector<double> b_n(ngroups * bs), b_het(ngroups * bs);

  for (size_t j0 = 0; j0 < m; j0 += bs) {
    // number of loci in this block
    size_t nb = std::min(bs, m - j0);

    // summarise each locus in the block by group
#pragma omp parallel num_threads(ncores)
{
    // per-thread buffers for the sums of each group
    std::vector<double> alt_sum(ngroups), n_alleles(ngroups), het2(ngroups);
#pragma omp for
    for (size_t jj = 0; jj < nb; jj++) {
      summ.locus(j0 + jj, alt_sum.data(), n_alleles.data(), het2.data());
      for (size_t g = 0; g < ngroups; g++) {
        size_t k = g * bs + jj;
        double freq_alt = alt_sum[g] / n_alleles[g];
        double freq_ref = 1 - freq_alt;
        b_freq_alt[k] = freq_alt;
        b_freq_ref[k] = freq_ref;
        b_n[k] = n_alleles[g];
        if (method_id == 0) {
          // for Hudson, we store the within population term
          // (freq_alt * freq_ref) / (n - 1) instead of het_obs (which is
          // not needed), so that it is computed once per group
          b_het[k] = (freq_alt * freq_ref) / (n_alleles[g] - 1);
        } else {
          b_het[k] = het2[g] / n_alleles[g];
        }
      }
    }
}

    // now compute fst for each pair of groups; each pair is handled by a
    // single thread
#pragma omp parallel num_threads(ncores)
{
    // per-thread buffers for the numerators and denominators of the block
    std::vector<double> num(bs), den(bs);
#pragma omp for schedule(dynamic, 16)
    for (size_t k = 0; k < n_combn; k++) {
      const double* fa1 = &b_freq_alt[pop1[k] * bs];
      const double* fa2 = &b_freq_alt[pop2[k] * bs];
      const double* fr1 = &b_freq_ref[pop1[k] * bs];
      const double* fr2 = &b_freq_ref[pop2[k] * bs];
      const double* n1 = &b_n[pop1[k] * bs];
      const double* n2 = &b_n[pop2[k] * bs];
      const double* h1 = &b_het[pop1[k] * bs];
      const double* h2 = &b_het[pop2[k] * bs];
      for (size_t jj = 0; jj < nb; jj++) {
        if (method_id == 0) {
          fst_hudson_num_den(fa1[jj], fr1[jj], h1[jj], fa2[jj], fr2[jj],
                             h2[jj], num[jj], den[jj]);
        } else if (method_id == 1) {
          fst_nei87_num_den(n1[jj], h1[jj], fa1[jj], fr1[jj],
                            n2[jj], h2[jj], fa2[jj], fr2[jj],
                            num[jj], den[jj]);
        } else {
          fst_wc84_num_den(n1[jj], fa1[jj], h1[jj], n2[jj], fa2[jj], h2[jj],
                           num[jj], den[jj]);
        }
      }
      consume(k, j0, nb, num.data(), den.data());
    }
}
  }
}

/******************************************************************************/
// Compute pairwise Fst for all pairs of groups, in a single pass over the FBM
// (see pairwise_fst_blocks()).
//
// method: one of "Hudson", "Nei87" or "WC84"
// pairwise_combn: a 2 x n_pairs matrix of (1-based) group indices
// by_locus: should Fst for each locus be returned?
// return_num_dem: should the numerator and denominator for each locus be
//   returned (instead of Fst by locus and the genome-wide estimate)?

// [[Rcpp::export]]
List pairwise_pop_fst_cpp(Environment BM,
                          const IntegerVector& rowInd,
                          const IntegerVector& colInd,
                          const IntegerVector& groupIds,
                          size_t ngroups,
                          const NumericVector& ploidy,
                          const IntegerMatrix& pairwise_combn,
                          std::string method,
                          bool by_locus,
                          bool return_num_dem,
                          int block_size,
                          int ncores) {

  // lookup tables and per-individual constants
  GroupedSummariesDipPseudo summ(BM, rowInd, colInd, groupIds, ngroups,
                                 ploidy);
  size_t m = summ.ncol(); // number of loci
  int method_id = fst_method_id(method);
  // pairs of groups (converted to 0-based indices)
  std::vector<int> pop1, pop2;
  fst_pairs(pairwise_combn, ngroups, pop1, pop2);
  size_t n_combn = pop1.size();

  // genome-wide sums of numerators and denominators for each pair
  std::vector<double> sum_num(n_combn, 0.0), sum_den(n_combn, 0.0);
  // if we need to return fst by locus, we create the matrix, else make an
  // empty one
  NumericMatrix fst_locus(by_locus ? m : 0, by_locus ? n_combn : 0);
  // if we need to return numerator and denominator, we use fst_locus for the
  // numerator and fst_locus_dem for the denominator
  NumericMatrix fst_locus_dem(return_num_dem ? m : 0,
                              return_num_dem ? n_combn : 0);
  double* p_fst_locus = REAL(fst_locus);
  double* p_fst_locus_dem = REAL(fst_locus_dem);

  pairwise_fst_blocks(summ, pop1, pop2, method_id, block_size, ncores,
    [&](size_t k, size_t j0, size_t nb, const double* num,
        const double* den) {
      // position of the first locus of the block in the output matrices
      size_t out_offset = k * m + j0;
      double s_num = sum_num[k], s_den = sum_den[k];
      for (size_t jj = 0; jj < nb; jj++) {
        if (by_locus) {
          if (!return_num_dem) {
            p_fst_locus[out_offset + jj] = num[jj] / den[jj];
          } else {
            p_fst_locus[out_offset + jj] = num[jj];
            p_fst_locus_dem[out_offset + jj] = den[jj];
          }
        }
        // only use loci for which both numerator and denominator are valid
        if (!std::isnan(num[jj]) && !std::isnan(den[jj])) {
          s_num += num[jj];
          s_den += den[jj];
        }
      }
      sum_num[k] = s_num;
      sum_den[k] = s_den;
    });

  if (!return_num_dem) {
    NumericVector fst_tot(n_combn);
    for (size_t k = 0; k < n_combn; k++) {
      fst_tot[k] = sum_num[k] / sum_den[k];
    }
    return List::create(Named("fst_locus") = fst_locus,
                        Named("fst_tot") = fst_tot);
  } else {
    return List::create(Named("Fst_by_locus_num") = fst_locus,
                        Named("Fst_by_locus_den") = fst_locus_dem);
  }
}

/******************************************************************************/
// Sums of the numerators and denominators of pairwise Fst in windows, for
// all pairs of groups, in a single pass over the FBM (see
// pairwise_fst_blocks()), without storing values by locus.
//
// For each locus (in the order of colInd), locus_first and locus_last give
// the first and last window (0-based) that include it (see
// .window_definitions(); a locus is in no window if locus_last <
// locus_first). For each window and pair, we accumulate the sum and the
// number of the non-missing numerators and denominators, in the order of the
// loci (as runner::mean_run() does when applied to each column of numerators
// and denominators). Returns a list of 4 matrices (windows x pairs): sum_num,
// n_num, sum_den and n_den.

// [[Rcpp::export]]
List windows_pairwise_pop_fst_cpp(Environment BM,
                                  const IntegerVector& rowInd,
                                  const IntegerVector& colInd,
                                  const IntegerVector& groupIds,
                                  size_t ngroups,
                                  const NumericVector& ploidy,
                                  const IntegerMatrix& pairwise_combn,
                                  std::string method,
                                  const IntegerVector& locus_first,
                                  const IntegerVector& locus_last,
                                  int n_windows,
                                  int block_size,
                                  int ncores) {

  // lookup tables and per-individual constants
  GroupedSummariesDipPseudo summ(BM, rowInd, colInd, groupIds, ngroups,
                                 ploidy);
  size_t m = summ.ncol(); // number of loci
  int method_id = fst_method_id(method);
  std::vector<int> pop1, pop2;
  fst_pairs(pairwise_combn, ngroups, pop1, pop2);
  size_t n_combn = pop1.size();
  if ((size_t) locus_first.size() != m || (size_t) locus_last.size() != m) {
    stop("locus_first and locus_last should have one value per locus");
  }
  for (size_t j = 0; j < m; j++) {
    if (locus_last[j] >= locus_first[j] &&
        (locus_first[j] < 0 || locus_last[j] >= n_windows)) {
      stop("invalid windows for a locus");
    }
  }
  const int* p_first = locus_first.begin();
  const int* p_last = locus_last.begin();

  // sums and counts for each window (rows) and pair (columns)
  NumericMatrix sum_num(n_windows, n_combn), n_num(n_windows, n_combn);
  NumericMatrix sum_den(n_windows, n_combn), n_den(n_windows, n_combn);
  double* p_sum_num = REAL(sum_num);
  double* p_n_num = REAL(n_num);
  double* p_sum_den = REAL(sum_den);
  double* p_n_den = REAL(n_den);

  pairwise_fst_blocks(summ, pop1, pop2, method_id, block_size, ncores,
    [&](size_t k, size_t j0, size_t nb, const double* num,
        const double* den) {
      // the windows of this pair
      size_t offset = k * (size_t) n_windows;
      for (size_t jj = 0; jj < nb; jj++) {
        size_t j = j0 + jj;
        for (int w = p_first[j]; w <= p_last[j]; w++) {
          if (!std::isnan(num[jj])) {
            p_sum_num[offset + w] += num[jj];
            p_n_num[offset + w] += 1;
          }
          if (!std::isnan(den[jj])) {
            p_sum_den[offset + w] += den[jj];
            p_n_den[offset + w] += 1;
          }
        }
      }
    });

  return List::create(_["sum_num"] = sum_num, _["n_num"] = n_num,
                      _["sum_den"] = sum_den, _["n_den"] = n_den);
}

/******************************************************************************/
