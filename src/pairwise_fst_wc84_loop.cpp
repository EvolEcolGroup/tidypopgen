#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List pairwise_fst_wc84_loop(NumericMatrix pairwise_combn,
                   NumericMatrix n,
                   NumericMatrix freq_alt,
                   NumericMatrix het_obs,
                   bool by_locus,
                   bool return_num_dem) {
  
  int n_loci = n.nrow();
  int n_combn = pairwise_combn.ncol();
  NumericVector fst_tot(n_combn); 
  
  // if we need to return fst by locus, we create the matrix, else make an empty one
  NumericMatrix fst_locus(by_locus ? n_loci : 0, by_locus ? n_combn : 0);
  // if we need to return numerator and denominator, we use fst_locus for the
  // numerator and fst_locus_dem for the denominator
  NumericMatrix fst_locus_dem(return_num_dem ? n_loci: 0, return_num_dem ? n_combn : 0);
  
  for (int i_col = 0; i_col < n_combn; ++i_col) {
    IntegerVector pops(pairwise_combn.nrow());
    for (int i = 0; i < pairwise_combn.nrow(); ++i) {
      pops[i] = pairwise_combn(i, i_col) - 1;  // Subtract 1 for 0-based indexing
    }
    int r = pops.size(); // usually 2
    
    NumericMatrix an(n_loci, r);
    NumericMatrix p(n_loci, r);
    NumericMatrix h(n_loci, r);
    
    for (int j = 0; j < r; ++j) {
      an(_, j) = n(_, pops[j]);
      p(_, j) = freq_alt(_, pops[j]);
      h(_, j) = het_obs(_, pops[j]);
    }
    
    NumericMatrix n_ind(n_loci, r);
    for (int i = 0; i < n_loci; ++i) {
      for (int j = 0; j < r; ++j) {
        n_ind(i, j) = an(i, j) / 2.0;
      }
    }
    
    NumericVector n_total(n_loci);
    NumericVector n_bar(n_loci);
    NumericVector n_c(n_loci);
    for (int i = 0; i < n_loci; ++i) {
      double sum_n = 0.0, sum_sq = 0.0;
      for (int j = 0; j < r; ++j) {
        sum_n += n_ind(i, j);
        sum_sq += pow(n_ind(i, j), 2);
      }
      n_total[i] = sum_n;
      n_bar[i] = sum_n / r;
      n_c[i] = (sum_n - sum_sq / sum_n) / (r - 1);
    }
    
    NumericVector p_bar(n_loci);
    NumericVector s2(n_loci);
    NumericVector h_bar(n_loci);
    
    for (int i = 0; i < n_loci; ++i) {
      double sum_pn = 0.0, sum_sq_diff = 0.0, sum_h = 0.0;
      for (int j = 0; j < r; ++j) {
        sum_pn += p(i, j) * n_ind(i, j);
        sum_h += h(i, j) * n_ind(i, j);
      }
      p_bar[i] = sum_pn / n_total[i];
      h_bar[i] = sum_h / n_total[i];
      for (int j = 0; j < r; ++j) {
        sum_sq_diff += pow(p(i, j) - p_bar[i], 2) * n_ind(i, j);
      }
      s2[i] = sum_sq_diff / (n_bar[i] * (r - 1));
    }
    
    NumericVector a(n_loci), b(n_loci), c(n_loci), numerator(n_loci), denominator(n_loci);
    for (int i = 0; i < n_loci; ++i) {
      a[i] = n_bar[i] / n_c[i] * (s2[i] -
        (1.0 / (n_bar[i] - 1.0)) *
        (p_bar[i] * (1 - p_bar[i]) - ((r - 1.0) / r) * s2[i] - h_bar[i] / 4.0));
      b[i] = n_bar[i] / (n_bar[i] - 1.0) *
        (p_bar[i] * (1 - p_bar[i]) - ((r - 1.0) / r) * s2[i] -
        ((2 * n_bar[i] - 1.0) / (4.0 * n_bar[i])) * h_bar[i]);
      c[i] = h_bar[i] / 2.0;
      
      numerator[i] = a[i];
      denominator[i] = a[i] + b[i] + c[i];
      
      if (by_locus) {
        if (!return_num_dem) {
          fst_locus(i, i_col) = numerator[i] / denominator[i];
        } else {
          fst_locus(i, i_col) = numerator[i];
          fst_locus_dem(i, i_col) = denominator[i];
        }
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
