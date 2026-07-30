/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// Function to count heterozygotes and na alleles in a matrix of genotypes
// by individual, generalised to arbitrary and mixed ploidy: a genotype is
// heterozygous whenever its dosage is neither 0 nor the individual's own
// ploidy (pseudohaploid dosages, which only ever take values 0 or 2, are
// correctly never counted as heterozygous under this same rule). Kept as a
// separate function from gt_ind_hetero() so the diploid fast path stays
// free of the extra ploidy lookup and comparison this requires.


// [[Rcpp::export]]
IntegerMatrix gt_ind_hetero_polyploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const NumericVector& ploidy,
                                   int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci

  IntegerMatrix het_counts(2, n); // matrix of 2 rows, n_het and n_na

  #pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x > -1){
        if (x > 0 && x < ploidy[i]){ // count heterozygote
          het_counts(0, i) += 1;
        }
      } else {
        // count missingness
        het_counts(1, i) += 1;
      }
    }
  }


  return het_counts;
}

/******************************************************************************/
