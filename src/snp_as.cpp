/******************************************************************************/

#include "pairwise_bitpacked.h"

/******************************************************************************/

// Increment the allele sharing counts for a block of loci.
// k: FBM (n x n) of tcrossprod(dos - 1), where dos is the dosage matrix with
//    missing values replaced by 1; i.e. the number of loci at which the two
//    individuals are the same homozygote, minus the number at which they are
//    opposite homozygotes
// k2: FBM (n x n) of tcrossprod(na), where na is 1 for valid genotypes and 0
//    for missing ones; i.e. the number of loci valid in both individuals
// byte_genotype: the genotype held by each raw byte (see
//    pairwise_byte_genotypes() in R)
// Genotypes are compared with the bit-packed kernel in pairwise_bitpacked.h

// [[Rcpp::export]]
void increment_as_counts(Environment k,
                         Environment k2,
                         Environment BM,
                         const IntegerVector& rowInd,
                         const IntegerVector& colInd,
                         const IntegerVector& byte_genotype,
                         int ncores) {

  std::size_t n = rowInd.size();
  double* K = pw_fbm_double_ptr(k, n);
  double* K2 = pw_fbm_double_ptr(k2, n);

  // counts for the packed lower triangle: (same homozygote, opposite
  // homozygotes, both valid)
  std::vector<std::size_t> tri_start;
  std::vector<uint32_t> counts = pw_bitpacked_counts<PwAlleleSharing>(
    BM, rowInd, colInd, byte_genotype, tri_start, ncores);

  // add them to both triangles of the (symmetric) accumulators
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = i; j < n; j++) {
      const uint32_t* c =
        &counts[(tri_start[i] + j - i) * PwAlleleSharing::N_COUNTS];
      double dos_cross = static_cast<double>(c[0]) - c[1];
      double n_valid = c[2];
      K[j + i * n] += dos_cross;
      K2[j + i * n] += n_valid;
      if (j != i) {
        K[i + j * n] += dos_cross;
        K2[i + j * n] += n_valid;
      }
    }
  }
}

/******************************************************************************/
