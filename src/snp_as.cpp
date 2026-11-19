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
// Genotypes are compared with the bit-packed kernel in pairwise_bitpacked.h

// [[Rcpp::export]]
void increment_as_counts(Environment k,
                         Environment k2,
                         Environment BM,
                         const IntegerVector& rowInd,
                         const IntegerVector& colInd,
                         int ncores) {

  std::size_t n = rowInd.size();
  double* K = pw_fbm_double_ptr(k, n);
  double* K2 = pw_fbm_double_ptr(k2, n);

  // counts for the packed lower triangle: (same homozygote, opposite
  // homozygotes, both valid)
  std::vector<std::size_t> tri_start;
  std::vector<uint32_t> counts = pw_bitpacked_counts<PwAlleleSharing>(
    BM, rowInd, colInd, tri_start, ncores);

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

// Counts for allele sharing within groups of individuals: rowInd lists the
// individuals ordered by group, and groupEnd gives, for each of them, the
// (1-based, inclusive) position in rowInd of the last individual of its
// group (as in king_counts_within_groups()). Only the pairs within each group
// are visited, packing the genotypes of all individuals once (rather than
// once for each group).
// Returns a matrix with a column for each pair (i, j), with i <= j in the
// same group (in the order of the packed triangle of each group, one group
// after the other), and 2 rows: the number of loci at which the two
// individuals are the same homozygote minus the number at which they are
// opposite homozygotes (the sum of the products of their centred dosages),
// and the number of loci valid in both. Allele sharing is
// 0.5 * (1 + row 1 / row 2).

// [[Rcpp::export]]
NumericMatrix as_counts_within_groups(Environment BM,
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
  std::vector<uint32_t> counts = pw_bitpacked_counts<PwAlleleSharing>(
    BM, rowInd, colInd, tri_start, ncores, row_end);
  std::size_t n_pairs = counts.size() / PwAlleleSharing::N_COUNTS;
  NumericMatrix res(2, n_pairs);
  for (std::size_t p = 0; p < n_pairs; p++) {
    const uint32_t* c = &counts[p * PwAlleleSharing::N_COUNTS];
    res(0, p) = static_cast<double>(c[0]) - c[1];
    res(1, p) = c[2];
  }
  return res;
}

/******************************************************************************/
