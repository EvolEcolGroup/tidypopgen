/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/
// Summarise a grouped gen_tibble
// this only works for diploid data


// [[Rcpp::export]]
ListOf<NumericMatrix> grouped_summaries_dip_pseudo_cpp(Environment BM,
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

  NumericMatrix freq(m, ngroups);
  NumericMatrix ref_freq(m, ngroups);
  NumericMatrix valid_alleles(m, ngroups);
  NumericMatrix heterozygotes(m, ngroups);
  
  // multiplier to convert dosages for pseudohaploid data
  NumericVector pseudo_multiplier (n);
  for (size_t i = 0; i < n; i++) {
    pseudo_multiplier[i] = 1 / (3 - ploidy[i]);
  }

#pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x > -1){
        freq(j, groupIds[i]) += x * pseudo_multiplier[i];
        valid_alleles(j, groupIds[i]) += ploidy[i];
        if (x==1){
          // we add 2 so that we can then divide by valid alleles
          // we don't have to worry about pseudohaploids, as heterozygote
          // counts are meaningless for them
          heterozygotes(j, groupIds[i]) +=2;
        }
      }
    }
    // now for each group, divide freq by valid_alleles
    for (size_t group_i = 0; group_i < ngroups; group_i++) {
      freq(j, group_i) = freq(j, group_i) / valid_alleles(j, group_i);
      ref_freq(j, group_i) = 1-freq(j, group_i);
      heterozygotes(j, group_i) = heterozygotes(j, group_i) / valid_alleles(j, group_i);
    }
  }

  return List::create(_["freq_alt"]  = freq,
                      _["freq_ref"]  = ref_freq,
                      _["n"] = valid_alleles,
                      _["het_obs"] = heterozygotes);
}

/******************************************************************************/
