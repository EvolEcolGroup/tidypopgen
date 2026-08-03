/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// Compute per-locus, per-group observed heterozygosity, defined simply as
// the proportion of heterozygous individuals (dosage neither 0 nor the
// individual's own ploidy) among the non-missing individuals in that group.
// Unlike the allele-frequency-based statistics (gene diversity, Fis, Fst,
// ...), this definition is ploidy-agnostic and requires no weighting
// decision, so it is valid even for groups that mix individuals of
// different ploidy.

// [[Rcpp::export]]
ListOf<NumericMatrix> grouped_het_obs_polyploid(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   const IntegerVector& groupIds,
                                   size_t ngroups,
                                   const NumericVector& ploidy,
                                   int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow(); // number of individuals
  size_t m = macc.ncol(); // number of loci

  NumericMatrix het_obs(m, ngroups);
  NumericMatrix n_valid(m, ngroups);

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x > -1){
        n_valid(j, groupIds[i]) += 1;
        if (x > 0 && x < ploidy[i]){
          het_obs(j, groupIds[i]) += 1;
        }
      }
    }
    for (size_t group_i = 0; group_i < ngroups; group_i++) {
      het_obs(j, group_i) = het_obs(j, group_i) / n_valid(j, group_i);
    }
  }

  return List::create(_["het_obs"] = het_obs,
                      _["n"] = n_valid);
}

/******************************************************************************/
