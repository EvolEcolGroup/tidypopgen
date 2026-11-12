#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

/******************************************************************************/
// Checks of the order of the loci table, used by is_loci_table_ordered() and
// gt_update_backingfile(). Each check is a single pass over the columns of the
// loci table, with a few values per chromosome rather than per locus, so that
// no column is copied or split (which, at tens of millions of loci, costs
// seconds and more than a GB of memory).
//
// Chromosomes are given as the factor `chromosome` of the loci table (its
// integer codes, with NA for missing chromosomes).

/******************************************************************************/
// Are values unsorted within any group of loci?
//
// This is equivalent to splitting x by group and checking whether any group
// is unsorted with is.unsorted(): a group of more than one locus with a
// missing value gives NA, so the result is TRUE if any group is unsorted, NA
// if no group is unsorted but some give NA, and FALSE otherwise.
//
// It also reports whether two consecutive non-missing values of a group are
// equal ("ties"). When all groups are sorted, this is the same as finding
// duplicated pairs of group and value.
//
// chromosome: the factor of chromosomes of the loci
// level_group: for each level of chromosome, the 1-based group of the loci
//   with that level (several levels can share a group); NA puts the loci in
//   the same group as the missing chromosomes
// Returns a named logical vector, c(unsorted, ties)

template <int RTYPE>
LogicalVector unsorted_by_group(const Vector<RTYPE>& x,
                                const IntegerVector& chromosome,
                                const IntegerVector& level_group) {
  R_xlen_t n_loci = x.size();
  if (chromosome.size() != n_loci) {
    stop("x and chromosome should have the same length");
  }
  int n_levels = level_group.size();
  // groups are 0 to n_levels - 1, and the missing chromosomes are group
  // n_levels
  int na_group = n_levels;
  // map each level (and NA, as level 0) to its 0-based group
  std::vector<int> group_of_code(n_levels + 1);
  group_of_code[0] = na_group;
  for (int l = 0; l < n_levels; l++) {
    int g = level_group[l];
    if (g == NA_INTEGER) {
      group_of_code[l + 1] = na_group;
    } else if (g < 1 || g > n_levels) {
      stop("level_group should be between 1 and the number of levels");
    } else {
      group_of_code[l + 1] = g - 1;
    }
  }

  // state of each group: number of loci, whether it has missing values,
  // whether it is unsorted, and its last non-missing value
  std::vector<R_xlen_t> n_in_group(n_levels + 1, 0);
  std::vector<char> has_na(n_levels + 1, 0);
  std::vector<char> unsorted(n_levels + 1, 0);
  std::vector<char> has_last(n_levels + 1, 0);
  std::vector<double> last(n_levels + 1, 0.0);
  bool ties = false;

  typedef typename traits::storage_type<RTYPE>::type stored_t;
  const stored_t* p_x = x.begin();
  const int* p_chr = chromosome.begin();
  for (R_xlen_t i = 0; i < n_loci; i++) {
    int code = p_chr[i];
    if (code == NA_INTEGER) {
      code = 0;
    } else if (code < 1 || code > n_levels) {
      stop("chromosome should be a factor with the levels of level_group");
    }
    int g = group_of_code[code];
    n_in_group[g]++;
    if (traits::is_na<RTYPE>(p_x[i])) {
      has_na[g] = 1;
      continue;
    }
    // compare with the previous value of the group (integers are exact as
    // doubles)
    double value = p_x[i];
    if (has_last[g]) {
      if (value < last[g]) {
        unsorted[g] = 1;
      } else if (value == last[g]) {
        ties = true;
      }
    }
    last[g] = value;
    has_last[g] = 1;
  }

  // combine the groups as any() over is.unsorted() would (is.unsorted() is
  // FALSE for a single value, even if it is missing)
  bool any_unsorted = false;
  bool any_na = false;
  for (int g = 0; g <= n_levels; g++) {
    if (n_in_group[g] < 2) continue;
    if (has_na[g]) {
      any_na = true;
    } else if (unsorted[g]) {
      any_unsorted = true;
    }
  }
  LogicalVector res = LogicalVector::create(
    _["unsorted"] = any_unsorted ? TRUE : (any_na ? NA_LOGICAL : FALSE),
    _["ties"] = ties
  );
  return res;
}

// [[Rcpp::export]]
LogicalVector unsorted_by_chromosome_cpp(SEXP x,
                                         const IntegerVector& chromosome,
                                         const IntegerVector& level_group) {
  switch (TYPEOF(x)) {
  case INTSXP:
    return unsorted_by_group<INTSXP>(IntegerVector(x), chromosome, level_group);
  case REALSXP:
    return unsorted_by_group<REALSXP>(NumericVector(x), chromosome,
                                      level_group);
  default:
    stop("x should be an integer or a numeric vector");
  }
}

/******************************************************************************/
// Is any chromosome split into more than one run of adjacent loci?
//
// This is equivalent to finding the runs of chromosome codes with rle() and
// checking whether any code appears in more than one run. As with rle(), each
// missing chromosome is a run of its own, so more than one missing chromosome
// gives TRUE.
//
// chromosome: the factor of chromosomes of the loci
// n_levels: the number of levels of chromosome

// [[Rcpp::export]]
bool chromosomes_not_adjacent_cpp(const IntegerVector& chromosome,
                                  int n_levels) {
  // whether each level has already been seen in an earlier run
  std::vector<char> seen(n_levels, 0);
  R_xlen_t n_missing = 0;
  // code of the previous locus (0 for none, as codes start at 1)
  int previous = 0;
  const int* p_chr = chromosome.begin();
  for (R_xlen_t i = 0; i < chromosome.size(); i++) {
    int code = p_chr[i];
    if (code == NA_INTEGER) {
      n_missing++;
      if (n_missing > 1) return true;
      previous = 0;
      continue;
    }
    if (code < 1 || code > n_levels) {
      stop("chromosome should be a factor with n_levels levels");
    }
    // a new run starts: its chromosome should not have been seen before
    if (code != previous) {
      if (seen[code - 1]) return true;
      seen[code - 1] = 1;
      previous = code;
    }
  }
  return false;
}

/******************************************************************************/
