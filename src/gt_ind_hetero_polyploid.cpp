/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <vector>

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

  // the outer loop is parallelised over loci, but the output is indexed by
  // individual, so different threads would otherwise race on the same
  // het_counts(*, i) cells; accumulate into thread-local counters (keeping
  // the loci-outer loop order, since bigstatsr's FBM is stored column-major
  // and this keeps access to macc(i, j) contiguous) and merge once per
  // thread rather than once per element.
  #pragma omp parallel num_threads(ncores)
  {
    std::vector<int> het_local(n, 0);
    std::vector<int> na_local(n, 0);

    #pragma omp for
    for (size_t j = 0; j < m; j++) {
      for (size_t i = 0; i < n; i++) {
        double x = macc(i, j);
        if (x > -1){
          if (x > 0 && x < ploidy[i]){ // count heterozygote
            het_local[i] += 1;
          }
        } else {
          // count missingness
          na_local[i] += 1;
        }
      }
    }

    #pragma omp critical
    {
      for (size_t i = 0; i < n; i++) {
        het_counts(0, i) += het_local[i];
        het_counts(1, i) += na_local[i];
      }
    }
  }

  return het_counts;
}

/******************************************************************************/
