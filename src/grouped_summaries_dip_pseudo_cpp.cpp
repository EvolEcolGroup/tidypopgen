/******************************************************************************/

#include "grouped_summaries_dip_pseudo.h"

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

  // lookup tables and per-individual constants (see
  // grouped_summaries_dip_pseudo.h)
  GroupedSummariesDipPseudo summ(BM, rowInd, colInd, groupIds, ngroups,
                                 ploidy);

  size_t m = summ.ncol(); // number of loci

  NumericMatrix freq(m, ngroups);
  NumericMatrix ref_freq(m, ngroups);
  NumericMatrix valid_alleles(m, ngroups);
  NumericMatrix heterozygotes(m, ngroups);
  // raw pointers to the outputs, so that we don't touch Rcpp objects within
  // the parallel region
  double* p_freq = REAL(freq);
  double* p_ref_freq = REAL(ref_freq);
  double* p_valid = REAL(valid_alleles);
  double* p_het = REAL(heterozygotes);

#pragma omp parallel num_threads(ncores)
{
  // per-thread buffers for the sums of each group
  std::vector<double> alt_sum(ngroups), n_alleles(ngroups), het2(ngroups);

#pragma omp for
  for (size_t j = 0; j < m; j++) {
    summ.locus(j, alt_sum.data(), n_alleles.data(), het2.data());
    // now for each group, divide freq by valid_alleles
    for (size_t group_i = 0; group_i < ngroups; group_i++) {
      size_t k = j + group_i * m;
      p_freq[k] = alt_sum[group_i] / n_alleles[group_i];
      p_ref_freq[k] = 1 - p_freq[k];
      p_valid[k] = n_alleles[group_i];
      p_het[k] = het2[group_i] / n_alleles[group_i];
    }
  }
}

  return List::create(_["freq_alt"]  = freq,
                      _["freq_ref"]  = ref_freq,
                      _["n"] = valid_alleles,
                      _["het_obs"] = heterozygotes);
}

/******************************************************************************/
