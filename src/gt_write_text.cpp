/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <fstream>
#include <string>
#include <vector>

/******************************************************************************/

// Helper to convert an element of a character vector into a std::string (NA
// is converted to "NA", as it would be by paste() in R)
inline std::string as_string_na(const CharacterVector& x, R_xlen_t i) {
  if (CharacterVector::is_na(x[i])) {
    return "NA";
  }
  return std::string(CHAR(STRING_ELT(x, i)));
}

// Helper to convert a character vector into a vector of std::string
inline std::vector<std::string> as_string_vector(const CharacterVector& x) {
  std::vector<std::string> res(x.size());
  for (R_xlen_t i = 0; i < x.size(); i++) {
    res[i] = as_string_na(x, i);
  }
  return res;
}

/******************************************************************************/

// Write the genotypes of a FBM.code256 to a text file, one line per locus
//
// The genotypes are taken from the raw bytes of the FBM, and each byte is
// written as the corresponding element of `lut` (a vector of 256 strings,
// built in R by decoding code256, so that the same function can write any
// text format). Each line starts with the matching element of `prefix` (if
// `prefix` is not empty), and the genotypes are separated by `sep`. Lines are
// appended to `file` (which is opened in text mode, as done by
// utils::write.table()).
//
// @param BM the FBM.code256
// @param rowInd the (1-based) rows of the FBM (individuals) to write
// @param colInd the (1-based) columns of the FBM (loci) to write
// @param lut a vector of 256 strings, one for each raw value
// @param prefix a vector of strings to start each line (one per locus), or an
//   empty vector
// @param sep the separator between genotypes
// @param file the file to append to
// [[Rcpp::export]]
void gt_write_text_by_locus(Environment BM,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd,
                            const CharacterVector& lut,
                            const CharacterVector& prefix,
                            const std::string& sep,
                            const std::string& file) {
  XPtr<FBM> xpBM = BM["address"];
  const unsigned char* fbm_ptr =
    static_cast<const unsigned char*>(xpBM->matrix());
  size_t n_fbm = xpBM->nrow(); // number of rows in the full FBM

  size_t n = rowInd.size(); // number of individuals to write
  size_t m = colInd.size(); // number of loci to write
  bool has_prefix = prefix.size() > 0;
  if (lut.size() != 256) {
    stop("'lut' should have 256 elements");
  }
  if (has_prefix && ((size_t) prefix.size() != m)) {
    stop("'prefix' should have one element per locus");
  }
  std::vector<std::string> lut_str = as_string_vector(lut);
  // convert the rows to 0-based indices once
  std::vector<size_t> rows(n);
  for (size_t i = 0; i < n; i++) {
    rows[i] = rowInd[i] - 1;
  }

  std::ofstream out(file.c_str(), std::ios::out | std::ios::app);
  if (!out.is_open()) {
    stop("could not open the file " + file);
  }
  std::string line;
  for (size_t j = 0; j < m; j++) {
    // pointer to the start of the column of this locus in the FBM
    const unsigned char* col_ptr = fbm_ptr + n_fbm * (colInd[j] - 1);
    line.clear();
    if (has_prefix) {
      line += as_string_na(prefix, j);
    }
    for (size_t i = 0; i < n; i++) {
      // the first genotype only needs a separator after a prefix
      if ((i > 0) || has_prefix) {
        line += sep;
      }
      line += lut_str[col_ptr[rows[i]]];
    }
    line += '\n';
    out.write(line.data(), line.size());
    // check for user interrupts every so often
    if ((j % 10000) == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  out.close();
  if (out.fail()) {
    stop("error when writing to the file " + file);
  }
}

/******************************************************************************/

// Write the genotypes of a FBM.code256 to a text file, one line per individual
//
// Each raw byte of the FBM is first mapped to a category with `byte_cat`
// (from 0 to n_cat - 1), and the genotype is written as the string in `lut`
// for that category. `lut` holds n_cat strings for each set of strings: if
// `locus_lut` is empty, all loci use the first set of strings; otherwise,
// `locus_lut` gives the (0-based) set of strings to use for each locus (this
// allows writing genotypes as alleles, e.g. "A G", with the alleles of each
// locus, without creating strings for every locus). Each line starts with the
// matching element of `prefix` (one per individual), and genotypes are
// separated by `sep`. Lines are appended to `file` (which is opened in text
// mode, as done by utils::write.table()).
//
// As the FBM is stored by column (locus), individuals are processed in blocks
// whose lines are built in memory (up to about `block_bytes` bytes), so that
// each block only needs a single pass through the FBM.
//
// @param BM the FBM.code256
// @param rowInd the (1-based) rows of the FBM (individuals) to write
// @param colInd the (1-based) columns of the FBM (loci) to write
// @param byte_cat a vector of 256 integers giving the category of each raw
//   value
// @param lut a vector of strings, with n_cat strings for each set
// @param n_cat the number of categories
// @param locus_lut the set of strings used by each locus, or an empty vector
// @param prefix a vector of strings to start each line (one per individual)
// @param sep the separator between genotypes (and after the prefix)
// @param file the file to append to
// @param block_bytes the approximate memory (in bytes) used to store lines
// [[Rcpp::export]]
void gt_write_text_by_indiv(Environment BM,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd,
                            const IntegerVector& byte_cat,
                            const CharacterVector& lut,
                            const int n_cat,
                            const IntegerVector& locus_lut,
                            const CharacterVector& prefix,
                            const std::string& sep,
                            const std::string& file,
                            const double block_bytes = 256e6) {
  XPtr<FBM> xpBM = BM["address"];
  const unsigned char* fbm_ptr =
    static_cast<const unsigned char*>(xpBM->matrix());
  size_t n_fbm = xpBM->nrow(); // number of rows in the full FBM

  size_t n = rowInd.size(); // number of individuals to write
  size_t m = colInd.size(); // number of loci to write
  bool by_locus = locus_lut.size() > 0;
  if (byte_cat.size() != 256) {
    stop("'byte_cat' should have 256 elements");
  }
  if (by_locus && ((size_t) locus_lut.size() != m)) {
    stop("'locus_lut' should have one element per locus");
  }
  if ((size_t) prefix.size() != n) {
    stop("'prefix' should have one element per individual");
  }
  std::vector<std::string> lut_str = as_string_vector(lut);
  std::vector<std::string> prefix_str = as_string_vector(prefix);
  size_t n_sets = lut_str.size() / n_cat;
  // check that all categories and sets of strings are valid
  for (int b = 0; b < 256; b++) {
    if ((byte_cat[b] < 0) || (byte_cat[b] >= n_cat)) {
      stop("'byte_cat' should be between 0 and n_cat - 1");
    }
  }
  for (size_t j = 0; by_locus && (j < m); j++) {
    if ((locus_lut[j] < 0) || ((size_t) locus_lut[j] >= n_sets)) {
      stop("'locus_lut' refers to a set of strings that does not exist");
    }
  }

  // estimate the maximum length of a line, to decide how many individuals
  // to process at once
  size_t max_width = 0;
  for (size_t k = 0; k < lut_str.size(); k++) {
    max_width = std::max(max_width, lut_str[k].size());
  }
  double line_bytes = (double) m * (sep.size() + max_width) + 1;
  size_t block_n = std::max(1.0, std::floor(block_bytes / line_bytes));
  block_n = std::min(block_n, std::max(n, (size_t) 1));

  std::ofstream out(file.c_str(), std::ios::out | std::ios::app);
  if (!out.is_open()) {
    stop("could not open the file " + file);
  }
  std::vector<std::string> lines(block_n);
  for (size_t first = 0; first < n; first += block_n) {
    size_t last = std::min(first + block_n, n);
    // start each line with the information on the individual
    for (size_t i = first; i < last; i++) {
      std::string& line = lines[i - first];
      line.clear();
      line.reserve(line_bytes + prefix_str[i].size());
      line += prefix_str[i];
    }
    // move through the loci, adding the genotypes of the individuals in
    // this block (which are contiguous in the column of the FBM)
    for (size_t j = 0; j < m; j++) {
      const unsigned char* col_ptr = fbm_ptr + n_fbm * (colInd[j] - 1);
      const std::string* this_lut =
        &lut_str[by_locus ? (size_t) n_cat * locus_lut[j] : 0];
      for (size_t i = first; i < last; i++) {
        std::string& line = lines[i - first];
        line += sep;
        line += this_lut[byte_cat[col_ptr[rowInd[i] - 1]]];
      }
    }
    // write the lines of this block
    for (size_t i = first; i < last; i++) {
      std::string& line = lines[i - first];
      line += '\n';
      out.write(line.data(), line.size());
    }
    Rcpp::checkUserInterrupt();
  }
  out.close();
  if (out.fail()) {
    stop("error when writing to the file " + file);
  }
}

/******************************************************************************/
