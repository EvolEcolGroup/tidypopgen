/******************************************************************************/

#include "fbm_byte_acc.h"

/******************************************************************************/

// Sums, in windows, of the quantities used to estimate Tajima's D for each
// group, in a single pass over the FBM, without storing pi or the number of
// valid allele copies by locus. Windows can overlap (windows_pop_tajimas_d()),
// be the blocks of a jackknife (jackknife_pop_tajimas_d()), or be a single
// window with all the loci (pop_tajimas_d()).
//
// Following Bailey et al. (2025), each locus x has its own sample size n_x,
// the number of valid allele copies of the group at that locus. The
// quantities summed are those of tajimas_d_locus_terms() (in
// R/pop_tajimas_d.R):
// - k_hat: pi at the locus (0 if n_x < 2);
// - seg: 1 if pi > 0 (a segregating site), 0 otherwise;
// - theta_w: 1 / a1(n_x) for a segregating site, 0 otherwise, where
//   a1(n) = sum_{i=1}^{n-1} 1/i;
// - n: n_x;
// - n_loci: 1.
// Loci with no valid allele copies (n_x = 0) are ignored. Pi and n_x are
// computed as in gt_grouped_pi_polyploid() (pseudohaploids contribute half
// their 0/2 dosage, and every individual as many allele copies as its
// ploidy), which for diploids is the same as gt_grouped_pi_diploid().
//
// For each locus j (in the order of colInd), locus_first[j] and
// locus_last[j] give the first and last window (0-based) that include it (see
// .window_definitions(); a locus is in no window if locus_last <
// locus_first). If both are empty, all loci are in window 0 (so that
// genome-wide sums do not need any vector with one value per locus).
//
// Each block of loci is summarised in parallel (keeping pi and n_x of the
// block for each group), and its loci are then added to their windows in the
// order of the loci, so the sums do not depend on the number of cores or on
// the block size. Sums are accumulated in long double, as R's sum() and
// colSums() do, so that they are identical to summing the values of
// tajimas_d_locus_terms() with colSums() (or window_sums_cpp() with
// long_double = TRUE).
//
// Returns a list of 5 matrices (windows x groups): k_hat, seg, theta_w, n and
// n_loci.

// [[Rcpp::export]]
List tajimas_d_window_sums_cpp(Environment BM,
                               const IntegerVector& rowInd,
                               const IntegerVector& colInd,
                               const IntegerVector& groupIds,
                               int ngroups,
                               const NumericVector& ploidy,
                               bool is_pseudohaploid,
                               const IntegerVector& locus_first,
                               const IntegerVector& locus_last,
                               int n_windows,
                               int ncores,
                               int block_size) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }
  std::vector<size_t> group_ids = check_group_ids(groupIds, n, ngroups);
  if (block_size < 1) stop("block_size should be a positive integer");
  if (n_windows < 1) stop("n_windows should be a positive integer");

  // the windows of each locus (or a single window for all loci)
  bool single_window = (locus_first.size() == 0 && locus_last.size() == 0);
  if (!single_window) {
    if ((size_t) locus_first.size() != m || (size_t) locus_last.size() != m) {
      stop("locus_first and locus_last should have one value per locus");
    }
    for (size_t j = 0; j < m; j++) {
      if (locus_last[j] >= locus_first[j] &&
          (locus_first[j] < 0 || locus_last[j] >= n_windows)) {
        stop("invalid windows for a locus");
      }
    }
  }
  const int* p_first = locus_first.begin();
  const int* p_last = locus_last.begin();

  // pseudohaploid individuals contribute half their (0/2) dosage, and every
  // individual contributes as many valid alleles as its ploidy
  std::vector<double> weight(n), copies(n);
  // the largest number of allele copies in a group
  std::vector<double> group_copies(ngroups, 0);
  for (size_t i = 0; i < n; i++) {
    weight[i] = (is_pseudohaploid && ploidy[i] == 1) ? 0.5 : 1;
    copies[i] = ploidy[i];
    group_copies[group_ids[i]] += copies[i];
  }
  size_t n_max = 0;
  for (size_t g = 0; g < (size_t) ngroups; g++) {
    n_max = std::max(n_max, (size_t) group_copies[g]);
  }
  // 1 / a1(n) for n = 0, ..., n_max (only used for n >= 2), with a1(n)
  // computed as in harmonic_sum() (i.e. R's cumsum(), which accumulates in
  // long double)
  std::vector<double> inv_a1(n_max + 1, 0.0);
  long double a1 = 0.0;
  for (size_t k = 2; k <= n_max; k++) {
    a1 += 1 / (double) (k - 1);
    inv_a1[k] = 1 / (double) a1;
  }

  // sums for each window and group (element [g * n_windows + w])
  size_t n_cells = (size_t) n_windows * ngroups;
  std::vector<long double> k_hat(n_cells, 0.0), seg(n_cells, 0.0),
    theta_w(n_cells, 0.0), n_sum(n_cells, 0.0), n_loci(n_cells, 0.0);

  // pi and number of valid allele copies of the loci of the current block,
  // for each group (element [(j - j0) * ngroups + g]; blocks start at
  // multiples of bs, so j - j0 = j % bs)
  size_t bs = std::min((size_t) block_size, std::max(m, (size_t) 1));
  std::vector<double> b_pi(bs * ngroups), b_valid(bs * ngroups);

  alt_valid_counts(acc, group_ids, ngroups, weight, copies, ncores, bs,
                   [&](size_t j, size_t g, double alt,
                       double this_valid_alleles) {
                     size_t k = (j % bs) * ngroups + g;
                     // pi is not defined with fewer than 2 valid allele
                     // copies, where the locus contributes 0 (as in pixy)
                     b_pi[k] = (this_valid_alleles > 1) ?
                       (alt * (this_valid_alleles - alt) /
                         (this_valid_alleles * (this_valid_alleles - 1) / 2)) :
                       0;
                     b_valid[k] = this_valid_alleles;
                   },
                   [&](size_t j0, size_t j1) {
                     // add the loci of the block to their windows, in order
                     for (size_t j = j0; j < j1; j++) {
                       int w_first = single_window ? 0 : p_first[j];
                       int w_last = single_window ? 0 : p_last[j];
                       for (size_t g = 0; g < (size_t) ngroups; g++) {
                         size_t k = (j - j0) * ngroups + g;
                         double n_x = b_valid[k];
                         // loci with no valid allele copies are ignored
                         if (!(n_x > 0)) continue;
                         double pi = b_pi[k];
                         double is_seg = (pi > 0) ? 1 : 0;
                         double theta = (pi > 0) ? inv_a1[(size_t) n_x] : 0;
                         for (int w = w_first; w <= w_last; w++) {
                           size_t cell = g * (size_t) n_windows + w;
                           k_hat[cell] += pi;
                           seg[cell] += is_seg;
                           theta_w[cell] += theta;
                           n_sum[cell] += n_x;
                           n_loci[cell] += 1;
                         }
                       }
                     }
                   });

  // copy the sums into matrices
  NumericMatrix k_hat_mat(n_windows, ngroups), seg_mat(n_windows, ngroups),
    theta_w_mat(n_windows, ngroups), n_mat(n_windows, ngroups),
    n_loci_mat(n_windows, ngroups);
  for (size_t cell = 0; cell < n_cells; cell++) {
    k_hat_mat[cell] = (double) k_hat[cell];
    seg_mat[cell] = (double) seg[cell];
    theta_w_mat[cell] = (double) theta_w[cell];
    n_mat[cell] = (double) n_sum[cell];
    n_loci_mat[cell] = (double) n_loci[cell];
  }
  return List::create(_["k_hat"] = k_hat_mat, _["seg"] = seg_mat,
                      _["theta_w"] = theta_w_mat, _["n"] = n_mat,
                      _["n_loci"] = n_loci_mat);
}

/******************************************************************************/
