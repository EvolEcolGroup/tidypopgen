/******************************************************************************/

#include "pairwise_bitpacked.h"

/******************************************************************************/

// Increment the KING-robust counts for a block of loci.
// k: FBM (n x n) of the KING numerator, i.e. the number of loci at which both
//    individuals are heterozygous, minus twice the number of loci at which
//    they are opposite homozygotes
// n_Aa_i: FBM (n x n) where element (i, j) is the number of loci at which i
//    is heterozygous and j has a valid genotype
// Genotypes are compared with the bit-packed kernel in pairwise_bitpacked.h

// [[Rcpp::export]]
void increment_king_numerator(Environment k,
                              Environment n_Aa_i,
                              Environment BM,
                              const IntegerVector& rowInd,
                              const IntegerVector& colInd,
                              int ncores) {

  std::size_t n = rowInd.size();
  double* K = pw_fbm_double_ptr(k, n);
  double* N_Aa_i = pw_fbm_double_ptr(n_Aa_i, n);

  // counts for the packed lower triangle: (het-het, opposite homozygotes,
  // i het & j valid, j het & i valid)
  std::vector<std::size_t> tri_start;
  std::vector<uint32_t> counts =
    pw_bitpacked_counts<PwKing>(BM, rowInd, colInd, tri_start, ncores);

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

// Counts for KING-robust within groups of individuals: rowInd lists the
// individuals ordered by group, and groupEnd gives, for each of them, the
// (1-based, inclusive) position in rowInd of the last individual of its
// group. Only the pairs within each group are visited, packing the genotypes
// of all individuals once (rather than once for each group).
// Returns a matrix with a column for each pair (i, j), with i <= j in the
// same group (in the order of the packed triangle of each group, one group
// after the other), and 4 rows: the number of loci where both are
// heterozygous, where they are opposite homozygotes, where i is
// heterozygous and j has a valid genotype, and where j is heterozygous and i
// has a valid genotype.

// [[Rcpp::export]]
NumericMatrix king_counts_within_groups(Environment BM,
                                        const IntegerVector& rowInd,
                                        const IntegerVector& colInd,
                                        const IntegerVector& groupEnd,
                                        int ncores) {
  std::size_t n = rowInd.size();
  if ((std::size_t) groupEnd.size() != n) {
    Rcpp::stop("groupEnd should have one value per individual.");
  }
  std::vector<std::size_t> row_end(n);
  for (std::size_t i = 0; i < n; i++) {
    row_end[i] = groupEnd[i]; // the exclusive 0-based end
  }
  std::vector<std::size_t> tri_start;
  std::vector<uint32_t> counts =
    pw_bitpacked_counts<PwKing>(BM, rowInd, colInd, tri_start, ncores,
                                row_end);
  // (a local copy, as the constructor takes a reference)
  const int n_counts = PwKing::N_COUNTS;
  std::size_t n_pairs = counts.size() / n_counts;
  NumericMatrix res(n_counts, n_pairs);
  for (std::size_t k = 0; k < counts.size(); k++) {
    res[k] = counts[k];
  }
  return res;
}

/******************************************************************************/
