/******************************************************************************/

#include "pairwise_bitpacked.h"
#include <cmath>

/******************************************************************************/

// Increment the IBS counts for a block of loci.
// k: FBM (n x n) of the counts of alleles identical by state, i.e. the sum
//    over loci valid in both individuals of 2 - |g_i - g_j|
// k2: FBM (n x n) of the number of valid alleles compared, i.e. 2 times the
//    number of loci valid in both individuals
// Genotypes are compared with the bit-packed kernel in pairwise_bitpacked.h

// [[Rcpp::export]]
void increment_ibs_counts(Environment k,
                          Environment k2,
                          Environment BM,
                          const IntegerVector& rowInd,
                          const IntegerVector& colInd,
                          int ncores) {

  std::size_t n = rowInd.size();
  double* K = pw_fbm_double_ptr(k, n);
  double* K2 = pw_fbm_double_ptr(k2, n);

  // counts for the packed lower triangle: (n both valid, sum |g_i - g_j|)
  std::vector<std::size_t> tri_start;
  std::vector<uint32_t> counts =
    pw_bitpacked_counts<PwIbs>(BM, rowInd, colInd, tri_start, ncores);

  // add them to both triangles of the (symmetric) accumulators
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = i; j < n; j++) {
      const uint32_t* c = &counts[(tri_start[i] + j - i) * PwIbs::N_COUNTS];
      double n_valid = c[0];
      double ibs = 2 * n_valid - c[1];
      K[j + i * n] += ibs;
      K2[j + i * n] += 2 * n_valid;
      if (j != i) {
        K[i + j * n] += ibs;
        K2[i + j * n] += 2 * n_valid;
      }
    }
  }
}

/******************************************************************************/

// As increment_ibs_counts(), for PLINK 1.9's weighted correction for missing
// genotypes (see snp_ibs()): also increments
// k3: FBM (n x n) of the sum of the weights of the loci missing in both
//    individuals (the diagonal holds the sum of the weights of the loci
//    missing in each individual)
// wts: the integer weights of the loci in colInd (stored as doubles, as they
//    can be larger than the largest R integer); their sum over all the loci
//    must be less than 2^32, so that the counts can not overflow

// [[Rcpp::export]]
void increment_ibs_counts_weighted(Environment k,
                                   Environment k2,
                                   Environment k3,
                                   Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const NumericVector& wts,
                                   int ncores) {

  std::size_t n = rowInd.size();
  double* K = pw_fbm_double_ptr(k, n);
  double* K2 = pw_fbm_double_ptr(k2, n);
  double* K3 = pw_fbm_double_ptr(k3, n);

  if (wts.size() != colInd.size()) {
    Rcpp::stop("wts should have one value per locus.");
  }
  // the weights as unsigned 32 bit integers
  std::vector<uint32_t> locus_wts(wts.size());
  double wts_sum = 0;
  for (R_xlen_t l = 0; l < wts.size(); l++) {
    if (!(wts[l] >= 0) || wts[l] != std::floor(wts[l])) {
      Rcpp::stop("wts should be non-negative integers.");
    }
    wts_sum += wts[l];
    locus_wts[l] = static_cast<uint32_t>(wts[l]);
  }
  if (wts_sum >= 4294967296.0) {
    Rcpp::stop("The sum of wts should be less than 2^32.");
  }

  // counts for the packed lower triangle: (n both valid, sum |g_i - g_j|,
  // sum of the weights of the loci missing in both)
  std::vector<std::size_t> tri_start;
  std::vector<uint32_t> counts = pw_bitpacked_counts<PwIbsWeighted>(
    BM, rowInd, colInd, tri_start, ncores, std::vector<std::size_t>(),
    locus_wts);

  // add them to both triangles of the (symmetric) accumulators
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = i; j < n; j++) {
      const uint32_t* c =
        &counts[(tri_start[i] + j - i) * PwIbsWeighted::N_COUNTS];
      double n_valid = c[0];
      double ibs = 2 * n_valid - c[1];
      double w_missing = c[2];
      K[j + i * n] += ibs;
      K2[j + i * n] += 2 * n_valid;
      K3[j + i * n] += w_missing;
      if (j != i) {
        K[i + j * n] += ibs;
        K2[i + j * n] += 2 * n_valid;
        K3[i + j * n] += w_missing;
      }
    }
  }
}

/******************************************************************************/
