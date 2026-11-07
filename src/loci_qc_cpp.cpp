/******************************************************************************/

#include "fbm_byte_acc.h"

/******************************************************************************/

// defined in hwe.cpp
double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2,
               uint32_t midp);
std::vector<uint64_t> genotype_count_lut(const FBMByteAcc& acc);

/******************************************************************************/

// Compute the statistics of qc_report_loci() (minor allele frequency,
// missingness and HWE exact test p-value) in a single pass over the
// genotypes, rather than reading them once for each statistic.
//
// For each locus, we count the genotypes with dosage 0, 1 and 2 and the
// missing genotypes, separately for each group and ploidy (ploidy one for
// pseudohaploid individuals, which store their allele as a 0/2 dosage), and
// derive:
// - maf: the minor allele frequency over all individuals, as loci_maf()
// - missingness: the proportion of missing genotypes, as loci_missingness()
// - hwe_p: the smallest HWE p-value across groups (computed only if hwe is
//   TRUE, otherwise NA), as loci_hwe()
// Since counts are exact, results are identical to those functions. This
// requires all non-missing codes to be 0, 1 or 2 (i.e. raw calls or imputed
// calls, but not dosages). Loci are processed in blocks of block_size,
// checking for user interrupts between blocks.

// [[Rcpp::export]]
List loci_qc_cpp(Environment BM,
                 const IntegerVector& rowInd,
                 const IntegerVector& colInd,
                 const IntegerVector& groupIds,
                 int ngroups,
                 const NumericVector& ploidy,
                 bool hwe,
                 uint32_t midp,
                 int ncores,
                 int block_size = 100000) {

  FBMByteAcc acc(BM, rowInd, colInd);

  size_t n = acc.nrow(); // number of individuals
  size_t m = acc.ncol(); // number of loci

  if (!acc.has_int_codes(2)) {
    stop("loci_qc_cpp requires genotypes coded as 0, 1, 2 or missing");
  }
  if (static_cast<size_t>(ploidy.size()) != n) {
    stop("ploidy must have one entry per individual (rowInd)");
  }
  std::vector<size_t> group_ids = check_group_ids(groupIds, n, ngroups);

  // distinct ploidy values (1 and/or 2), and the class of each individual
  std::vector<double> ploidy_values;
  std::vector<size_t> ploidy_ids(n);
  for (size_t i = 0; i < n; i++) {
    size_t c = std::find(ploidy_values.begin(), ploidy_values.end(),
                         ploidy[i]) - ploidy_values.begin();
    if (c == ploidy_values.size()) ploidy_values.push_back(ploidy[i]);
    ploidy_ids[i] = c;
  }
  size_t nclasses = ploidy_values.size();
  // multiplier to convert dosages for pseudohaploid data (as in
  // alt_freq_dip_pseudo_cpp())
  std::vector<double> pseudo_multiplier(nclasses);
  for (size_t c = 0; c < nclasses; c++) {
    pseudo_multiplier[c] = 1 / (3 - ploidy_values[c]);
  }

  // genotypes are counted for each combination of group and ploidy class
  std::vector<size_t> count_ids(n);
  for (size_t i = 0; i < n; i++) {
    count_ids[i] = group_ids[i] * nclasses + ploidy_ids[i];
  }
  std::vector<uint64_t> lut = genotype_count_lut(acc);
  GroupedRows grouped(acc.rows(), count_ids, std::vector<size_t>(n, 0),
                      max_run_length<4>(lut));

  NumericVector maf(m), missingness(m), hwe_p(m, NA_REAL);
  double* p_maf = maf.begin();
  double* p_missingness = missingness.begin();
  double* p_hwe = hwe_p.begin();

  lut_sums<4>(acc, grouped, lut, ngroups * nclasses, ncores, block_size,
              [&](size_t j, const uint64_t* tot) {
                double n_alt = 0, n_valid = 0;
                uint64_t n_na = 0;
                double min_p = R_PosInf;
                for (size_t g = 0; g < static_cast<size_t>(ngroups); g++) {
                  uint64_t n0 = 0, n1 = 0, n2 = 0;
                  for (size_t c = 0; c < nclasses; c++) {
                    const uint64_t* cnt = tot + 4 * (g * nclasses + c);
                    // alternate alleles and valid alleles (sums of
                    // integers and halves, so they are exact)
                    n_alt += (cnt[1] + 2 * cnt[2]) * pseudo_multiplier[c];
                    n_valid += (cnt[0] + cnt[1] + cnt[2]) * ploidy_values[c];
                    n_na += cnt[3];
                    n0 += cnt[0];
                    n1 += cnt[1];
                    n2 += cnt[2];
                  }
                  if (hwe) {
                    min_p = std::min(min_p, SNPHWE2(n1, n0, n2, midp));
                  }
                }
                // minor allele frequency (NA if no valid alleles)
                double freq = (n_valid > 0) ? n_alt / n_valid : NA_REAL;
                if (!ISNAN(freq) && freq > 0.5) freq = 1 - freq;
                p_maf[j] = freq;
                p_missingness[j] = static_cast<double>(n_na) / n;
                if (hwe) p_hwe[j] = min_p;
              });

  return List::create(_["maf"] = maf,
                      _["missingness"] = missingness,
                      _["hwe_p"] = hwe_p);
}

/******************************************************************************/
