/******************************************************************************/

#include "gt_ind_counts.h"

/******************************************************************************/

// Function to count heterozygotes and na alleles in a matrix of genotypes
// by individual
// this only works for diploid data; see gt_ind_hetero_polyploid() for the
// generalised, ploidy-aware version used for polyploid/mixed-ploidy data.
// Both share the same lookup-table kernel (gt_ind_het_na_counts(), see
// gt_ind_counts.h); here all individuals use a single table, built from
// code256 with the same rules as decoding each genotype (a genotype is
// missing unless its decoded value is > -1, and heterozygous if it equals 1).


// [[Rcpp::export]]
IntegerMatrix gt_ind_hetero(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   int ncores) {

  // build the lookup table from the raw byte to (het, na) increments
  NumericVector code256 = BM["code256"];
  std::vector<uint64_t> lut(256);
  for (int b = 0; b < 256; b++) {
    double x = code256[b];
    bool valid = x > -1; // false for NA (as NaN)
    lut[b] = gt_pack_het_na(valid && (x == 1), !valid);
  }
  // a single table for all individuals
  std::vector<std::size_t> row_class(rowInd.size(), 0);

  return gt_ind_het_na_counts(BM, rowInd, colInd, lut, row_class, ncores);
}

/******************************************************************************/
