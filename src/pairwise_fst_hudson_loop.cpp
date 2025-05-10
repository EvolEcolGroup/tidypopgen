#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List pairwise_fst_hudson_loop(NumericMatrix pairwise_combn, List pop_freqs_df, bool by_locus, bool return_num_dem) {
  int ncol_combn = pairwise_combn.ncol();
  int n_loci = as<NumericMatrix>(pop_freqs_df["freq_alt"]).nrow();
  // variable used to define size of fst_locus, dummy value of 1 if we don't use by locus
  int n_loci_cols = 1;
  if (by_locus){
    n_loci_cols = n_loci;
  }
  NumericMatrix fst_locus(n_loci_cols, ncol_combn);
  // if we need to return numerator and denominator, we use fst_locus for the
  // numerator and fst_locus_dem for the denominator
  NumericMatrix fst_locus_dem(n_loci_cols, ncol_combn);
  NumericVector fst_tot(ncol_combn);

  NumericMatrix freq_alt = as<NumericMatrix>(pop_freqs_df["freq_alt"]);
  NumericMatrix freq_ref = as<NumericMatrix>(pop_freqs_df["freq_ref"]);
  NumericMatrix n = as<NumericMatrix>(pop_freqs_df["n"]);

  for (int i_col = 0; i_col < ncol_combn; ++i_col) {
    int pop1 = pairwise_combn(0, i_col) - 1; // Convert from R 1-based index to C++ 0-based index
    int pop2 = pairwise_combn(1, i_col) - 1;

    NumericVector numerator = pow(freq_alt(_, pop1) - freq_alt(_, pop2), 2) -
      (freq_alt(_, pop1) * freq_ref(_, pop1)) / (n(_, pop1) - 1) -
      (freq_alt(_, pop2) * freq_ref(_, pop2)) / (n(_, pop2) - 1);

    NumericVector denominator = freq_alt(_, pop1) * freq_ref(_, pop2) +
      freq_alt(_, pop2) * freq_ref(_, pop1);

    if (by_locus) {
      if (!return_num_dem) {
        fst_locus(_, i_col) = numerator / denominator;
      } else {
        fst_locus(_, i_col) = numerator;
        fst_locus_dem(_, i_col) = denominator;
      }
    }
    
    double mean_num = 0.0, mean_denom = 0.0;
    int count = 0;
    for (int i = 0; i < n_loci; ++i) {
      if (!NumericVector::is_na(numerator[i]) && !NumericVector::is_na(denominator[i])) {
        mean_num += numerator[i];
        mean_denom += denominator[i];
        count++;
      }
    }
    fst_tot[i_col] = mean_num / mean_denom;
  }

  if (!return_num_dem) {
    return List::create(Named("fst_locus") = fst_locus,
                        Named("fst_tot") = fst_tot);
  } else {
    return List::create(Named("Fst_by_locus_num") = fst_locus,
                        Named("Fst_by_locus_den") = fst_locus_dem);
  }

}
