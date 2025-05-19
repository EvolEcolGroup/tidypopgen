#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List pairwise_fst_nei87_loop(NumericMatrix pairwise_combn,
                      NumericMatrix n,
                      NumericMatrix het_obs,
                      NumericMatrix freq_alt,
                      NumericMatrix freq_ref,
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
    int pop1 = pairwise_combn(0, i_col) -1; // Convert from R 1-based index to C++ 0-based index
    int pop2 = pairwise_combn(1, i_col) -1;

    NumericMatrix n_pair(n_loci, 2), sHo(n_loci, 2), freqA(n_loci, 2), freqR(n_loci, 2);

    for (int i = 0; i < n_loci; ++i) {
      n_pair(i, 0) = n(i, pop1) / 2.0;
      n_pair(i, 1) = n(i, pop2) / 2.0;

      sHo(i, 0) = het_obs(i, pop1);
      sHo(i, 1) = het_obs(i, pop2);

      freqA(i, 0) = freq_alt(i, pop1);
      freqA(i, 1) = freq_alt(i, pop2);
      freqR(i, 0) = freq_ref(i, pop1);
      freqR(i, 1) = freq_ref(i, pop2);
    }

    NumericVector mHo(n_loci), np(n_loci), mn(n_loci), msp2(n_loci), mp2(n_loci), mHs(n_loci), Ht(n_loci);
    NumericVector Dst(n_loci), Dstp(n_loci), Htp(n_loci), fst(n_loci);

    for (int i = 0; i < n_loci; ++i) {
      // Sample size and het_obs means
      int valid = 0;
      double nsum = 0.0, inv_nsum = 0.0;
      double ho_sum = 0.0;

      for (int j = 0; j < 2; ++j) {
        if (!NumericVector::is_na(n_pair(i, j))) {
          valid++;
          double nij = n_pair(i, j);
          ho_sum += sHo(i, j);
          nsum += 1.0;
          inv_nsum += 1.0 / nij;
        }
      }

      np[i] = valid;
      mn[i] = (inv_nsum > 0) ? nsum / inv_nsum : NA_REAL;
      mHo[i] = ho_sum / valid;

      double sp2a = pow(freqA(i, 0), 2) + pow(freqA(i, 1), 2);
      double sp2r = pow(freqR(i, 0), 2) + pow(freqR(i, 1), 2);
      double sp2 = sp2a + sp2r;
      msp2[i] = sp2 / 2.0;
      double freqAm = (freqA(i, 0) + freqA(i, 1)) / 2.0;
      double freqRm = (freqR(i, 0) + freqR(i, 1)) / 2.0;
      mp2[i] = pow(freqAm, 2) + pow(freqRm, 2);

//      double hs1 = 1.0 - sp2 - sHo(i, 0) / (2.0 * n_pair(i, 0));
//      double hs2 = 1.0 - sp2 - sHo(i, 1) / (2.0 * n_pair(i, 1));

      mHs[i] = mn[i] / (mn[i] - 1.0) * (1.0 - msp2[i] - mHo[i] / (2.0 * mn[i]));
      Ht[i] = 1.0 - mp2[i] + mHs[i] / (mn[i] * np[i]) - mHo[i] / (2.0 * mn[i] * np[i]);
      Dst[i] = Ht[i] - mHs[i];
      Dstp[i] = (np[i] / (np[i] - 1.0)) * Dst[i];
      Htp[i] = mHs[i] + Dstp[i];

      if (by_locus) {
        if (!return_num_dem) {
          fst_locus(i, i_col) =  Dstp[i] / Htp[i];
        } else {
          fst_locus(i, i_col) = Dstp[i];
          fst_locus_dem(i, i_col) = Htp[i];
        }
      }


    }

    double mean_num = 0.0, mean_denom = 0.0;
    int count = 0;
    for (int i = 0; i < n_loci; ++i) {
      if (!NumericVector::is_na(Dstp[i]) && !NumericVector::is_na(Htp[i])) {
        mean_num += Dstp[i];
        mean_denom += Htp[i];
        count++;
      }
    }
    fst_tot[i_col] = mean_num / mean_denom;

    // fst_tot[i_col] = mean(Dstp) / mean(Htp);
  }

  if (!return_num_dem) {
    return List::create(Named("fst_locus") = fst_locus,
                        Named("fst_tot") = fst_tot);
  } else {
    return List::create(Named("Fst_by_locus_num") = fst_locus,
                        Named("Fst_by_locus_den") = fst_locus_dem);
  }
}
