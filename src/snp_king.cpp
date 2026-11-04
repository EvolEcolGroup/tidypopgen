/******************************************************************************/

#include "pairwise_bitpacked.h"

/******************************************************************************/

// Increment the KING-robust counts for a block of loci.
// k: FBM (n x n) of the KING numerator, i.e. the number of loci at which both
//    individuals are heterozygous, minus twice the number of loci at which
//    they are opposite homozygotes
// n_Aa_i: FBM (n x n) where element (i, j) is the number of loci at which i
//    is heterozygous and j has a valid genotype
// byte_genotype: the genotype held by each raw byte (see
//    pairwise_byte_genotypes() in R)
// Genotypes are compared with the bit-packed kernel in pairwise_bitpacked.h

// [[Rcpp::export]]
void increment_king_numerator(Environment k,
                              Environment n_Aa_i,
                              Environment BM,
                              const IntegerVector& rowInd,
                              const IntegerVector& colInd,
                              const IntegerVector& byte_genotype,
                              int ncores) {

  std::size_t n = rowInd.size();
  double* K = pw_fbm_double_ptr(k, n);
  double* N_Aa_i = pw_fbm_double_ptr(n_Aa_i, n);

  // counts for the packed lower triangle: (het-het, opposite homozygotes,
  // i het & j valid, j het & i valid)
  std::vector<std::size_t> tri_start;
  std::vector<uint32_t> counts =
    pw_bitpacked_counts<PwKing>(BM, rowInd, colInd, byte_genotype, tri_start, ncores);

  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = i; j < n; j++) {
      const uint32_t* c = &counts[(tri_start[i] + j - i) * PwKing::N_COUNTS];
      double numerator = static_cast<double>(c[0]) - 2.0 * c[1];
      // the numerator is symmetric
      K[j + i * n] += numerator;
      if (j != i) {
        K[i + j * n] += numerator;
      }
      // the counts of heterozygotes are not: (i, j) and (j, i)
      N_Aa_i[i + j * n] += c[2];
      if (j != i) {
        N_Aa_i[j + i * n] += c[3];
      }
    }
  }
}

/******************************************************************************/
