/******************************************************************************/

#include "pairwise_bitpacked.h"

/******************************************************************************/

// Increment the IBS counts for a block of loci.
// k: FBM (n x n) of the counts of alleles identical by state, i.e. the sum
//    over loci valid in both individuals of 2 - |g_i - g_j|
// k2: FBM (n x n) of the number of valid alleles compared, i.e. 2 times the
//    number of loci valid in both individuals
// byte_genotype: the genotype held by each raw byte (see
//    pairwise_byte_genotypes() in R)
// Genotypes are compared with the bit-packed kernel in pairwise_bitpacked.h

// [[Rcpp::export]]
void increment_ibs_counts(Environment k,
                          Environment k2,
                          Environment BM,
                          const IntegerVector& rowInd,
                          const IntegerVector& colInd,
                          const IntegerVector& byte_genotype,
                          int ncores) {

  std::size_t n = rowInd.size();
  double* K = pw_fbm_double_ptr(k, n);
  double* K2 = pw_fbm_double_ptr(k2, n);

  // counts for the packed lower triangle: (n both valid, sum |g_i - g_j|)
  std::vector<std::size_t> tri_start;
  std::vector<uint32_t> counts =
    pw_bitpacked_counts<PwIbs>(BM, rowInd, colInd, byte_genotype, tri_start, ncores);

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
