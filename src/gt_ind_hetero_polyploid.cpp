/******************************************************************************/

#include "gt_ind_counts.h"
#include <map>

/******************************************************************************/

// Function to count heterozygotes and na alleles in a matrix of genotypes
// by individual, generalised to arbitrary and mixed ploidy: a genotype is
// heterozygous whenever its dosage is neither 0 nor the individual's own
// ploidy (pseudohaploid dosages, which only ever take values 0 or 2, are
// correctly never counted as heterozygous under this same rule).
// It shares the lookup-table kernel of gt_ind_hetero() (see
// gt_ind_counts.h), with one table per distinct ploidy value, so the
// ploidy comparison is done once per (ploidy, byte) rather than once per
// genotype.


// [[Rcpp::export]]
IntegerMatrix gt_ind_hetero_polyploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const NumericVector& ploidy,
                                   int ncores) {

  if (ploidy.size() != rowInd.size()) {
    Rcpp::stop("ploidy should have one value per individual");
  }
  NumericVector code256 = BM["code256"];

  // assign each distinct ploidy value to a lookup table ("class"), and
  // record for each individual the offset of its table
  std::map<double, std::size_t> class_of_ploidy;
  std::vector<std::size_t> row_class(rowInd.size());
  std::vector<uint64_t> lut;
  for (R_xlen_t i = 0; i < ploidy.size(); i++) {
    double this_ploidy = ploidy[i];
    auto found = class_of_ploidy.find(this_ploidy);
    if (found == class_of_ploidy.end()) {
      // new ploidy: append its table of 256 entries
      std::size_t offset = lut.size();
      class_of_ploidy[this_ploidy] = offset;
      for (int b = 0; b < 256; b++) {
        double x = code256[b];
        bool valid = x > -1; // false for NA (as NaN)
        lut.push_back(
          gt_pack_het_na(valid && (x > 0 && x < this_ploidy), !valid));
      }
      row_class[i] = offset;
    } else {
      row_class[i] = found->second;
    }
  }

  return gt_ind_het_na_counts(BM, rowInd, colInd, lut, row_class, ncores);
}

/******************************************************************************/
