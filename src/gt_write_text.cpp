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
