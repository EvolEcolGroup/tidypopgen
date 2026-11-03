/******************************************************************************/

#include "fbm_byte_acc.h"

/******************************************************************************/

// Shared by alt_freq_dip_pseudo_cpp() and alt_freq_polyploid_cpp(): if
// as_counts = TRUE, return a matrix with two columns, the count of alternate
// alleles and the count of valid alleles; otherwise, a vector of frequencies,
// with na_freq for loci with no valid alleles.
static SEXP alt_freq_output(const FBMByteAcc& acc,
                            const std::vector<double>& weight,
                            const std::vector<double>& copies,
                            bool as_counts,
                            double na_freq,
                            int ncores,
                            int block_size) {
  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci
  // all individuals belong to a single group
  std::vector<size_t> group_ids(n, 0);

  if (as_counts) {
    // matrix with the counts of alternate alleles and of valid alleles
    NumericMatrix freq_mat(m, 2);
    double* n_alt = freq_mat.begin();
    double* n_valid = n_alt + m;
    alt_valid_counts(acc, group_ids, 1, weight, copies, ncores, block_size,
                     [&](size_t j, size_t g, double alt, double valid) {
                       n_alt[j] = alt;
                       n_valid[j] = valid;
                     });
    colnames(freq_mat) = CharacterVector::create("n_alt", "n_valid");
    return freq_mat;
  }

  // only return the frequencies, to avoid allocating the counts
  NumericVector freq(m);
  double* p_freq = freq.begin();
  alt_valid_counts(acc, group_ids, 1, weight, copies, ncores, block_size,
                   [&](size_t j, size_t g, double alt, double valid) {
                     p_freq[j] = (valid > 0) ? alt / valid : na_freq;
                   });
  return freq;
}

/******************************************************************************/

// Estimate the alternate allele frequency at each locus for diploid and
// pseudohaploid data. If as_counts = TRUE, it returns a matrix with two
// columns, the count of alternate alleles and the count of valid alleles;
// otherwise, a vector of frequencies (NA when there are no valid alleles).
// Loci are processed in blocks of block_size, checking for user interrupts
// between blocks.

// [[Rcpp::export]]
SEXP alt_freq_dip_pseudo_cpp(Environment BM,
                             const IntegerVector& rowInd,
                             const IntegerVector& colInd,
                             const NumericVector& ploidy,
                             int ncores,
                             bool as_counts,
                             int block_size = 100000) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals

  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }

  // multiplier to convert dosages for pseudohaploid data, and number of
  // valid alleles contributed by each individual
  std::vector<double> pseudo_multiplier(n), copies(n);
  for (size_t i = 0; i < n; i++) {
    pseudo_multiplier[i] = 1 / (3 - ploidy[i]);
    copies[i] = ploidy[i];
  }

  return alt_freq_output(acc, pseudo_multiplier, copies, as_counts, NA_REAL,
                         ncores, block_size);
}

/******************************************************************************/

// Estimate the alternate allele frequency at each locus for polyploid (and
// mixed ploidy) data, where each individual contributes its dosage and as
// many valid alleles as its ploidy. Frequencies are NaN when there are no
// valid alleles.

// [[Rcpp::export]]
SEXP alt_freq_polyploid_cpp(Environment BM,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd,
                            const NumericVector& ploidy,
                            int ncores,
                            bool as_counts,
                            int block_size = 100000) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals

  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }

  // dosages are counted as they are, and each individual contributes as many
  // valid alleles as its ploidy
  std::vector<double> weight(n, 1), copies(n);
  for (size_t i = 0; i < n; i++) {
    copies[i] = ploidy[i];
  }

  return alt_freq_output(acc, weight, copies, as_counts, R_NaN, ncores,
                         block_size);
}

/******************************************************************************/
