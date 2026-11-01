#ifndef TIDYPOPGEN_GROUPED_SUMMARIES_DIP_PSEUDO_H
#define TIDYPOPGEN_GROUPED_SUMMARIES_DIP_PSEUDO_H

/******************************************************************************/

#include <bigstatsr/BMAcc.h>
#include <algorithm>
#include <vector>

/******************************************************************************/
// Summarise the genotypes of diploid (and pseudohaploid) individuals by group,
// one locus at a time.
//
// This is shared by pop_stats_dip_pseudo_cpp() and pairwise_pop_fst_cpp(),
// which process loci in blocks. Rather than decoding each genotype through the
// code256 vector into a double (and accumulating directly into large
// loci x groups matrices), we read the raw bytes of the FBM and use small
// lookup tables (one entry per possible byte) to get the dosage, whether the
// genotype is valid, and whether it is a heterozygote. Sums are accumulated
// into per-group buffers (which stay in cache), in the same order (by
// individual) as in the original implementation, so results are identical.
//
// All Rcpp objects are converted to plain C++ containers in the constructor,
// so that locus() can be safely called from multiple OpenMP threads.

class GroupedSummariesDipPseudo {
public:
  GroupedSummariesDipPseudo(Environment BM,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd,
                            const IntegerVector& groupIds,
                            size_t ngroups,
                            const NumericVector& ploidy) {
    XPtr<FBM> xpBM = BM["address"];
    _pMat = static_cast<const unsigned char*>(xpBM->matrix());
    _nrow_fbm = xpBM->nrow();
    // convert the 1-based R indices to 0-based indices (checking the bounds)
    _row_ind = vec_int_to_size(rowInd, xpBM->nrow(), 1);
    _col_ind = vec_int_to_size(colInd, xpBM->ncol(), 1);
    _ngroups = ngroups;

    size_t n = _row_ind.size();
    if (groupIds.size() != (R_xlen_t) n || ploidy.size() != (R_xlen_t) n) {
      stop("groupIds and ploidy should have one value per individual");
    }
    _group.resize(n);
    _ploidy.resize(n);
    _multiplier.resize(n);
    for (size_t i = 0; i < n; i++) {
      if (groupIds[i] < 0 || (size_t) groupIds[i] >= ngroups) {
        stop("groupIds should be between 0 and ngroups - 1");
      }
      _group[i] = groupIds[i];
      _ploidy[i] = ploidy[i];
      // multiplier to convert dosages for pseudohaploid data
      _multiplier[i] = 1 / (3 - ploidy[i]);
    }

    // lookup tables for each possible value of the raw byte
    NumericVector code256 = BM["code256"];
    _dosage.resize(256);
    _is_valid.resize(256);
    _het_inc.resize(256);
    for (int b = 0; b < 256; b++) {
      double x = code256[b];
      _dosage[b] = x;
      // NA (and NaN) fail this comparison, so they are not valid
      _is_valid[b] = (x > -1);
      // we add 2 for heterozygotes so that we can then divide by valid
      // alleles (heterozygote counts are meaningless for pseudohaploids)
      _het_inc[b] = (x == 1) ? 2.0 : 0.0;
    }
  }

  size_t nrow() const { return _row_ind.size(); }
  size_t ncol() const { return _col_ind.size(); }
  size_t ngroups() const { return _ngroups; }

  // Summarise locus j (0-based index into colInd). For each group g, on return
  // alt_sum[g] is the sum of alternative allele dosages (corrected for
  // pseudohaploids), n_alleles[g] the number of valid alleles, and het2[g]
  // twice the number of heterozygotes. Each buffer must have ngroups elements.
  inline void locus(size_t j, double* alt_sum, double* n_alleles,
                    double* het2) const {
    std::fill(alt_sum, alt_sum + _ngroups, 0.0);
    std::fill(n_alleles, n_alleles + _ngroups, 0.0);
    std::fill(het2, het2 + _ngroups, 0.0);
    // pointer to the column of the FBM for this locus
    const unsigned char* col = _pMat + _col_ind[j] * _nrow_fbm;
    size_t n = _row_ind.size();
    for (size_t i = 0; i < n; i++) {
      unsigned char b = col[_row_ind[i]];
      if (_is_valid[b]) {
        int g = _group[i];
        alt_sum[g] += _dosage[b] * _multiplier[i];
        n_alleles[g] += _ploidy[i];
        // adding 0 for homozygotes leaves the sum unchanged
        het2[g] += _het_inc[b];
      }
    }
  }

private:
  const unsigned char* _pMat;
  size_t _nrow_fbm;
  std::vector<size_t> _row_ind;
  std::vector<size_t> _col_ind;
  size_t _ngroups;
  std::vector<int> _group;
  std::vector<double> _ploidy;
  std::vector<double> _multiplier;
  std::vector<double> _dosage;
  std::vector<char> _is_valid;
  std::vector<double> _het_inc;
};

/******************************************************************************/

#endif // TIDYPOPGEN_GROUPED_SUMMARIES_DIP_PSEUDO_H
