/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <vector>

/******************************************************************************/

// Function to count heterozygotes and na alleles in a matrix of genotypes
// by individual
// this only works for diploid data; see gt_ind_hetero_polyploid() for the
// generalised, ploidy-aware version used for polyploid/mixed-ploidy data.
// Kept separate (rather than folded into a single ploidy-aware function) so
// that the diploid hot path stays branch- and lookup-free.


// [[Rcpp::export]]
IntegerMatrix gt_ind_hetero(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
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
          if (x == 1){ // count heterozygote
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
