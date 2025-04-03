#include <Rcpp.h>
using namespace Rcpp;

// Function to compute mean without using sugar to deal with na removal
double mean_cpp(NumericVector x) {
  int n = x.size();
  double sum = 0.0;
  int count = 0;
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(x[i])) {
      sum += x[i];
      count++;
    }
  }
  return count > 0 ? sum / count : NA_REAL;
}

// [[Rcpp::export]]
List pairwise_fst_hudson_loop(NumericMatrix pairwise_combn, List pop_freqs_df, bool by_locus) {
  int ncol_combn = pairwise_combn.ncol();
  int n_loci = as<NumericMatrix>(pop_freqs_df["freq_alt"]).nrow();
  int n_loci_cols; // variable used to define size of fst_locus
  if (by_locus){
    n_loci_cols = n_loci;
  } else {
    n_loci_cols = 1; // if not by locus, we only need one dummy column to avoid assigning a lot of memory
  }
  NumericMatrix fst_locus(n_loci_cols, ncol_combn);
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
      for (int i = 0; i < n_loci; ++i) {
        fst_locus(i, i_col) = numerator[i] / denominator[i];
      }
    }
    
    fst_tot[i_col] = mean_cpp(numerator) / mean_cpp(denominator);
  }
  
  return List::create(Named("fst_locus") = fst_locus,
                      Named("fst_tot") = fst_tot);
}