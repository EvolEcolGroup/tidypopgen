/******************************************************************************/

#include "fbm_byte_acc.h"

/******************************************************************************/

// Compute per-locus, per-group observed heterozygosity, defined simply as
// the proportion of heterozygous individuals (dosage neither 0 nor the
// individual's own ploidy) among the non-missing individuals in that group.
// Unlike the allele-frequency-based statistics (gene diversity, Fis, Fst,
// ...), this definition is ploidy-agnostic and requires no weighting
// decision, so it is valid even for groups that mix individuals of
// different ploidy. Loci are processed in blocks of block_size, checking for
// user interrupts between blocks.

// [[Rcpp::export]]
ListOf<NumericMatrix> grouped_het_obs_polyploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   size_t ngroups,
                                   const NumericVector& ploidy,
                                   int ncores,
                                   int block_size = 100000) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }
  std::vector<size_t> group_ids = check_group_ids(groupIds, n, ngroups);

  // one lookup table per distinct ploidy: the low 32 bits count valid
  // genotypes and the high 32 bits heterozygous ones
  std::vector<double> ploidy_values;
  std::vector<size_t> lut_ids(n);
  for (size_t i = 0; i < n; i++) {
    size_t c = std::find(ploidy_values.begin(), ploidy_values.end(),
                         ploidy[i]) - ploidy_values.begin();
    if (c == ploidy_values.size()) ploidy_values.push_back(ploidy[i]);
    lut_ids[i] = c;
  }
  std::vector<uint64_t> luts(256 * ploidy_values.size(), 0);
  for (size_t c = 0; c < ploidy_values.size(); c++) {
    for (int b = 0; b < 256; b++) {
      double x = acc.code(b);
      if (x > -1) {
        uint64_t het = (x > 0 && x < ploidy_values[c]) ? 1 : 0;
        luts[256 * c + b] = 1 | (het << 32);
      }
    }
  }
  GroupedRows grouped(acc.rows(), group_ids, lut_ids,
                      max_run_length<2>(luts));

  NumericMatrix het_obs(m, ngroups);
  NumericMatrix n_valid(m, ngroups);
  double* p_het = het_obs.begin();
  double* p_valid = n_valid.begin();

  lut_sums<2>(acc, grouped, luts, ngroups, ncores, block_size,
              [&](size_t j, const uint64_t* tot) {
                for (size_t g = 0; g < ngroups; g++) {
                  double this_valid = static_cast<double>(tot[2 * g]);
                  p_het[j + g * m] =
                    static_cast<double>(tot[2 * g + 1]) / this_valid;
                  p_valid[j + g * m] = this_valid;
                }
              });

  return List::create(_["het_obs"] = het_obs,
                      _["n"] = n_valid);
}

/******************************************************************************/
