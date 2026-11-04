#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <unordered_map>
#include <zlib.h>
#include <bigstatsr/BMCodeAcc.h>

using namespace Rcpp;
using namespace std;

// Function to read a line from a gzipped file
inline bool gz_getline(gzFile file, std::string &line) {
  char buffer[1024];
  line.clear();
  while (gzgets(file, buffer, sizeof(buffer)) != Z_NULL) {
    line.append(buffer);
    if (line.back() == '\n') {
      line.pop_back();
      return true;
    }
  }
  return !line.empty();
}


// Function to remove the carriage return at the end of a line (CRLF endings)
inline void strip_carriage_return(std::string &line) {
  if (!line.empty() && line.back() == '\r') {
    line.pop_back();
  }
}

// Function to check whether a file name ends in .gz (ignoring case)
inline bool is_gz_file(const std::string &filename) {
  if (filename.size() < 3) {
    return false;
  }
  std::string ext = filename.substr(filename.size() - 3);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  return ext == ".gz";
}

// Function to move a pointer from the start of a field to the start of the
// field n_fields further along a tab separated line (it returns NULL if the
// line has fewer fields)
inline const char* skip_fields(const char* p, size_t n_fields) {
  for (size_t i = 0; (i < n_fields) && (p != NULL); i++) {
    p = std::strchr(p, '\t');
    if (p != NULL) {
      p++;
    }
  }
  return p;
}

// Function to get the ploidy of a genotype (the number of alleles before the
// first ':' or tab), without splitting the line. A missing genotype "." gives
// no information on ploidy, so it returns 0
inline int genotype_ploidy(const char* p) {
  if ((p[0] == '.') && ((p[1] == ':') || (p[1] == '\t') || (p[1] == '\0'))) {
    return 0;
  }
  int n_alleles = 1;
  for (; (*p != ':') && (*p != '\t') && (*p != '\0'); p++) {
    // count the genotype delimiters (unphased or phased)
    if ((*p == '/') || (*p == '|')) {
      n_alleles++;
    }
  }
  return n_alleles;
}


// Function to convert a vector of characters into an R character vector,
// reusing a single R string for each distinct character (rather than creating
// a string for each element)
inline CharacterVector chars_to_strings(const std::vector<char> &x) {
  CharacterVector res(x.size());
  // the R string of each of the 256 possible characters (created when first
  // needed; each is protected by being stored in res)
  std::vector<SEXP> char_strings(256, R_NilValue);
  for (size_t i = 0; i < x.size(); i++) {
    unsigned char c = x[i];
    if (char_strings[c] == R_NilValue) {
      char_strings[c] = Rf_mkCharLen(&x[i], 1);
    }
    SET_STRING_ELT(res, i, char_strings[c]);
  }
  return res;
}

// Function to split a string by multiple delimiters and return a vector of tokens
inline std::vector<std::string> split_line(const std::string &s, const std::string &delimiters) {
  std::vector<std::string> tokens;
  size_t start = s.find_first_not_of(delimiters), end = 0;

  while ((end = s.find_first_of(delimiters, start)) != std::string::npos) {
    tokens.push_back(s.substr(start, end - start));
    start = s.find_first_not_of(delimiters, end);
  }
  if (start != std::string::npos) {
    tokens.push_back(s.substr(start));
  }
  return tokens;
}

// Function to split the metadata for a locus and return a vector of tokens
// this only needs to read the first 5 elements in the string
inline std::vector<std::string> split_locus_meta(const std::string &s,
                                                 const std::string &delimiters) {
  std::vector<std::string> tokens;
  size_t start = s.find_first_not_of(delimiters), end = 0;

  for (int i = 0; i < 5; i++) {
    end = s.find_first_of(delimiters, start);
    tokens.push_back(s.substr(start, end - start));
    start = s.find_first_not_of(delimiters, end);
  }
  return tokens;
}



/******************************************************************************
 * Create a loci table and get ploidy and number of individuals from vcf
 *
 * This function parses the VCF for biallelic snps, and creates a loci table
 * with information for those snps. It also creates a boolean vector to indicate
 * which loci are biallelic, and a vector with the ploidy of each individual.
 *
 * To limit memory use with many loci, chromosomes are returned as an integer
 * index (1-based) into a vector of chromosome names (in order of appearance),
 * alleles (a single character for biallelic snps) are stored as characters,
 * and locus names are only read if they are needed.
 *
 * @param filename The name of the VCF file to parse
 * @param read_names Whether to read the names (IDs) of the loci
 */
// [[Rcpp::export]]
List vcf_loci_table(std::string filename, bool read_names = true) {
  // open the file
  bool is_gz = is_gz_file(filename);
  std::ifstream vcfFile;
  gzFile vcfGzFile;

  if (is_gz) {
    vcfGzFile = gzopen(filename.c_str(), "rb");
    if (vcfGzFile == NULL) {
      stop("Error opening gzipped file");
    }
  } else {
    vcfFile.open(filename);
    if (!vcfFile.is_open()) {
      stop("Error opening file");
    }
  }

  // is a locus biallelic
  std::vector<bool> biallelic;
  biallelic.reserve(10e6);
  // loci metadata (preallocate 10 million loci, as a start)
  // chromosomes are stored as an index into the vector of their names (with a
  // lookup table to find the index of a name)
  std::vector<int> chromosome_index;
  chromosome_index.reserve(10e6);
  std::vector<std::string> chromosome_names;
  std::unordered_map<std::string, int> chromosome_lookup;
  int this_chromosome = -1; // index of the chromosome of the previous locus
  std::vector<std::string> marker_id;
  if (read_names) {
    marker_id.reserve(10e6);
  }
  std::vector<int32_t> physical_pos;
  physical_pos.reserve(10e6);
  // alleles of biallelic snps are a single character
  std::vector<char> allele1;
  allele1.reserve(10e6);
  std::vector<char> allele2;
  allele2.reserve(10e6);
  // individual metadata
  std::vector<int> ploidy;
  std::vector<std::string> sample_names;


  std::string line;
  std::string line_old; // used to get the samples from the last line of the header
  int n_header_lines =0;

  bool get_sample_names = true;
  bool get_ploidy = true;
  // individuals whose genotypes have all been missing (".") so far, and the
  // number of variants in which we still look for their ploidy (we only look
  // in the first few variants, so that the search does not slow down
  // parsing)
  std::vector<size_t> unknown_ploidy;
  int ploidy_lines_left = 1000;

  // Read the file line by line until the end of file
  while (true) {
    bool gotLine;
    if (is_gz) {
      gotLine = gz_getline(vcfGzFile, line);
    } else {
      gotLine = (bool)std::getline(vcfFile, line);
    }

    // Check if we reached the end of the file
    if (!gotLine) {
      break;
    }
    // remove the carriage return of CRLF line endings, and skip empty lines
    strip_carriage_return(line);
    if (line.empty()) {
      continue;
    }

    // Skip header lines
    if (line[0] == '#') {
      line_old = line; // save the last line of the header
      n_header_lines++;
      continue;
    }

    // Check if we need to get the sample names (this is only the case for
    // the first line after headers)
    if (get_sample_names) {
      std::vector<std::string> fields = split_line(line_old, "\t");
      for (unsigned int i = 9; i < fields.size(); ++i) {
        sample_names.push_back(fields[i]);
      }
      get_sample_names = false;
    }

    // Extract chromosome, position, and
    // Split the VCF line into fields (only for the first few fields)
    std::vector<std::string> fields = split_locus_meta(line, "\t");

    // Check if the marker is biallelic SNP (REF and ALT allele are single character)
    // if yes, get the metadata for this locus
    if ((fields[3].length() ==1) && (fields[4].length() ==1)) {
      biallelic.push_back(true);
      // find the index of the chromosome (loci are usually sorted by
      // chromosome, so we only need a lookup when the chromosome changes)
      if ((this_chromosome < 0) ||
          (fields[0] != chromosome_names[this_chromosome])) {
        std::unordered_map<std::string, int>::iterator found =
          chromosome_lookup.find(fields[0]);
        if (found == chromosome_lookup.end()) {
          this_chromosome = chromosome_names.size();
          chromosome_lookup[fields[0]] = this_chromosome;
          chromosome_names.push_back(fields[0]);
        } else {
          this_chromosome = found->second;
        }
      }
      chromosome_index.push_back(this_chromosome + 1);
      physical_pos.push_back(std::stoi(fields[1]));
      if (read_names) {
        marker_id.push_back(fields[2]);
      }
      allele1.push_back(fields[4][0]);
      allele2.push_back(fields[3][0]);
    } else {
      biallelic.push_back(false);
    }


    // Check if we need to get ploidy (this is the case for the first line
    // of real data; a missing genotype "." gives no information on ploidy, so
    // for those individuals we use their first genotype that is not ".")
    if (get_ploidy) {
      // the genotypes start at the 10th field
      const char* p = skip_fields(line.c_str(), 9);
      while (p != NULL) {
        int this_ploidy = genotype_ploidy(p);
        if (this_ploidy == 0) {
          unknown_ploidy.push_back(ploidy.size());
        }
        ploidy.push_back(this_ploidy);
        p = skip_fields(p, 1);
      }
      get_ploidy = false;
    } else if (!unknown_ploidy.empty() && (ploidy_lines_left > 0)) {
      ploidy_lines_left--;
      const char* p = line.c_str();
      size_t field = 0; // the field that p points to
      size_t n_unknown = 0; // individuals still with unknown ploidy
      for (size_t k = 0; k < unknown_ploidy.size(); k++) {
        size_t i = unknown_ploidy[k];
        // move to the genotype of this individual (the field 9 + i)
        if (p != NULL) {
          p = skip_fields(p, i + 9 - field);
          field = i + 9;
        }
        int this_ploidy = (p != NULL) ? genotype_ploidy(p) : 0;
        if (this_ploidy > 0) {
          ploidy[i] = this_ploidy;
        } else {
          unknown_ploidy[n_unknown++] = i;
        }
      }
      unknown_ploidy.resize(n_unknown);
    }

  }

  // individuals with only missing genotypes (in the variants where we looked
  // for their ploidy) get the largest ploidy of the other individuals (or 2 if
  // no ploidy could be estimated)
  if (!unknown_ploidy.empty()) {
    int default_ploidy = 0;
    for (size_t i = 0; i < ploidy.size(); ++i) {
      default_ploidy = std::max(default_ploidy, ploidy[i]);
    }
    if (default_ploidy == 0) {
      default_ploidy = 2;
    }
    for (size_t i = 0; i < ploidy.size(); ++i) {
      if (ploidy[i] == 0) {
        ploidy[i] = default_ploidy;
      }
    }
  }

  // close the file
  if (is_gz) {
    gzclose(vcfGzFile);
  } else {
    vcfFile.close();
  }

  // create the loci table; alleles are converted to strings reusing a single
  // R string for each character, and the C++ vectors are freed as soon as
  // they have been copied into R
  IntegerVector chromosome_r(chromosome_index.begin(), chromosome_index.end());
  std::vector<int>().swap(chromosome_index);
  IntegerVector physical_pos_r(physical_pos.begin(), physical_pos.end());
  std::vector<int32_t>().swap(physical_pos);
  CharacterVector allele1_r = chars_to_strings(allele1);
  std::vector<char>().swap(allele1);
  CharacterVector allele2_r = chars_to_strings(allele2);
  std::vector<char>().swap(allele2);
  List loci_table = List::create(Named("chromosome") = chromosome_r,
                                 _["physical.pos"] = physical_pos_r,
                                 _["allele1"] = allele1_r,
                                 _["allele2"] = allele2_r);
  if (read_names) {
    loci_table["marker.ID"] = marker_id;
    std::vector<std::string>().swap(marker_id);
  }

  // Return a list containing the loci table and the ploidy vector
  return (List::create(Named("loci_tbl") = loci_table,
                       Named("chromosome_names") = chromosome_names,
                       Named("sample_names") = sample_names,
                       Named("ploidy") = ploidy,
                       Named("biallelic") = biallelic,
                       Named("n_header_lines") = n_header_lines));
}


// Function to count the number of alternate alleles from genotype information
//
/******************************************************************************
 * count the number of alternate alleles from genotype information
 *
 * missingValue is max_ploidy +1
 *
 * @param filename The name of the VCF file to parse
 * This function parses a line of the VCF and counts the number of alternate alleles,
 * storing it in an array of raw values (than can then be written directly to the
 * file backed matrix).
 *
 * @param string line The line of the VCF file to parse
 * @param unsigned char* alt_counts The array to store the counts
 * @param int n The number of individuals
 * @param int missingValue The value to use for missing genotypes
 * @return false if a genotype has more alleles than max_ploidy
 */

// find the position of the field separators of a line of the VCF (there are
// 8 + n of them, and separator_pos[8 + n] is set to the end of the line);
// returns false if the line has fewer than 8 + n separators (i.e. fewer
// genotypes than individuals)
inline bool find_separators(const std::string &line,
                            size_t* separator_pos,
                            size_t n,
                            const std::string &separator) {
  bool complete = true;
  separator_pos[0] = line.find(separator, 0);
  complete &= (separator_pos[0] != std::string::npos);
  for (size_t i = 1; i < 8+n; i++) {
    separator_pos[i] = line.find(separator, separator_pos[i-1]+1);
    complete &= (separator_pos[i] != std::string::npos);
  }
  separator_pos[8+n] = line.length();
  return complete;
}

inline bool count_alt_alleles(const std::string &line,
                              unsigned char* alt_counts,
                              size_t* separator_pos,
                              size_t n,
                              const size_t missingValue,
                              const std::string &separator
) {
  // find the position of all separators
  find_separators(line, separator_pos, n, separator);

  // now move through each of the genotypes and process them
  // we ignore the first 9 fields (so start from position 9)
  const char* line_ptr = line.c_str();
  for (size_t i=0; i<n; i++){
    // the genotypes start at pos +1
    const char* p = line_ptr + separator_pos[i+8] + 1;
    const char* field_end = line_ptr + separator_pos[i+9];
    // for a biallelic locus, each allele is a single character ('0', '1' or
    // '.'), and alleles are delimited by '/' or '|' (unphased or phased), so
    // we jump from one allele to the next, until we reach the end of the
    // genotype section (':' or the end of the field)
    size_t n_alt = 0;
    size_t n_alleles = 0;
    bool is_missing = false;
    while (true) {
      n_alt += (p[0] == '1');
      // if any allele is missing, the genotype is missing
      is_missing |= (p[0] == '.');
      n_alleles++;
      if ((p + 1 >= field_end) || ((p[1] != '/') && (p[1] != '|'))) {
        break;
      }
      p += 2;
    }
    // a genotype can not have more alleles than max_ploidy
    if (n_alleles > missingValue - 1) {
      return false;
    }
    alt_counts[i] = is_missing ? missingValue : n_alt;
  }
  return true;
}

/******************************************************************************
 * phased bytes from the genotypes of a line of the VCF
 *
 * For each individual, this stores the byte that encodes its phased diploid
 * genotype (see `?code_phased`): phased_lut[2 * a + b] for a genotype a|b
 * (with a and b equal to 0 for the reference allele and 1 for the alternate
 * allele), or missing_value if any of its alleles is missing ('.'), or if it
 * is a single missing allele ('.'). The phase of a homozygote is known, so
 * homozygotes can also be unphased ('0/0' and '1/1').
 *
 * @return -1 if all genotypes could be read; otherwise the (0-based) index of
 *   the first individual whose genotype could not be stored as a phased
 *   genotype, with `error` set to 1 if the genotype is not diploid, or to 2 if
 *   it is an unphased heterozygote ('0/1' or '1/0'); or 0, with `error` set
 *   to 3, if the line has fewer genotypes than individuals.
 */
inline long phased_genotype_bytes(const std::string &line,
                                  unsigned char* bytes,
                                  size_t* separator_pos,
                                  size_t n,
                                  const unsigned char missing_value,
                                  const unsigned char* phased_lut,
                                  const std::string &separator,
                                  int &error) {
  if (!find_separators(line, separator_pos, n, separator)) {
    // the line is too short: there are fewer genotypes than individuals
    error = 3;
    return 0;
  }
  const char* line_ptr = line.c_str();
  for (size_t i = 0; i < n; i++) {
    const char* p = line_ptr + separator_pos[i+8] + 1;
    const char* field_end = line_ptr + separator_pos[i+9];
    // the end of the genotype (the first ':' or the end of the field)
    const char* gt_end = p;
    while ((gt_end < field_end) && (*gt_end != ':')) {
      gt_end++;
    }
    size_t len = gt_end - p;
    if ((len == 1) && (p[0] == '.')) {
      // a single missing allele
      bytes[i] = missing_value;
    } else if ((len == 3) && ((p[1] == '|') || (p[1] == '/'))) {
      // a diploid genotype (each allele of a biallelic locus is a single
      // character)
      if ((p[0] == '.') || (p[2] == '.')) {
        bytes[i] = missing_value;
      } else {
        int a = (p[0] == '1');
        int b = (p[2] == '1');
        if ((p[1] == '/') && (a != b)) {
          error = 2;
          return i;
        }
        bytes[i] = phased_lut[2 * a + b];
      }
    } else {
      error = 1;
      return i;
    }
  }
  return -1;
}

/**
 * @brief Reads genotype data from a VCF file and stores alternate allele counts in a file-backed matrix.
 *
 * For each biallelic locus in the VCF file, counts the number of alternate alleles per individual and writes these counts to the provided file-backed matrix (FBM). Supports both plain text and gzipped VCF files. Skips the specified number of header lines and processes only loci marked as biallelic. Verifies that no genotype exceeds the maximum ploidy inferred from the first variant line.
 *
 * @param filename Path to the VCF file (plain or gzipped).
 * @param biallelic Logical vector indicating which loci are biallelic and should be processed.
 * @param missing_value Value to assign for missing genotypes.
 * @param n_header_lines Number of header lines to skip in the VCF file.
 * @param phased_lut An empty vector to store alternate allele counts, or, to
 * store phased diploid genotypes, the bytes of the genotypes 0|0, 0|1, 1|0 and
 * 1|1 (in this order; see `phased_genotype_bytes()`).
 * @param sample_names The names of the individuals (used in error messages).
 * @return true if all genotypes are successfully read and stored.
 *
 * @throws Rcpp::exception if the file cannot be opened, if the end of file is reached prematurely, or if a genotype exceeds the maximum allowed ploidy (or, for phased genotypes, if a genotype is not diploid or is an unphased heterozygote).
 */

// [[Rcpp::export]]
bool vcf_genotypes_to_fbm(std::string filename,
                          Environment BM,
                          LogicalVector& biallelic,
                          const int missing_value,
                          const int n_header_lines,
                          const IntegerVector& phased_lut,
                          const CharacterVector& sample_names) {
  // are we storing phased genotypes?
  const bool phased = (phased_lut.size() > 0);
  unsigned char phased_lut_bytes[4] = {0, 0, 0, 0};
  if (phased) {
    if (phased_lut.size() != 4) {
      stop("'phased_lut' should be empty or have 4 elements");
    }
    for (int k = 0; k < 4; k++) {
      if ((phased_lut[k] < 0) || (phased_lut[k] > 255)) {
        stop("'phased_lut' should contain bytes (0 to 255)");
      }
      phased_lut_bytes[k] = (unsigned char) phased_lut[k];
    }
    if ((missing_value < 0) || (missing_value > 255)) {
      stop("'missing_value' should be a byte (0 to 255)");
    }
  }
  // open the file
  bool is_gz = is_gz_file(filename);
  std::ifstream vcfFile;
  gzFile vcfGzFile;

  if (is_gz) {
    vcfGzFile = gzopen(filename.c_str(), "rb");
    if (vcfGzFile == NULL) {
      stop("Error opening gzipped file");
    }
  } else {
    vcfFile.open(filename);
    if (!vcfFile.is_open()) {
      stop("Error opening file");
    }
  }

  // set up access to the file backed matrix
  XPtr<FBM_RW> xpBM = BM["address_rw"];
  unsigned char* ptr = static_cast<unsigned char*>(xpBM->matrix());
  //  const unsigned char* code_ptr;
  int n = xpBM->nrow(); // number of individuals


  // total number of loci in the vcf
  int n_loci = biallelic.size(); // number of loci

  // string to store the line read from the file
  std::string line;
  std::string tab_delimiter = "\t";

  // Skip the header
  for (int i = 0; i < n_header_lines; i++) {
    if (is_gz) {
      gz_getline(vcfGzFile, line);
    } else {
      std::getline(vcfFile, line);
    }
  }
  // array to store the genotypes
//  unsigned char* alt_counts = new unsigned char[n];
  std::unique_ptr<unsigned char[]> alt_counts(new unsigned char[n]);
  // position of the separator, the last value is the end of the string
//  size_t* separator_pos = new size_t[n+9];
  std::unique_ptr<size_t[]> separator_pos(new size_t[n+9]);

  // Read the genotypes line by line until the end of file
  for (int i_geno = 0; i_geno < n_loci; i_geno++) {
    bool gotLine;
    // skip empty lines (and any header line left), as in vcf_loci_table
    do {
      if (is_gz) {
        gotLine = gz_getline(vcfGzFile, line);
      } else {
        gotLine = (bool)std::getline(vcfFile, line);
      }
      if (gotLine) {
        strip_carriage_return(line);
      }
    } while (gotLine && (line.empty() || (line[0] == '#')));
    // Check if we reached the end of the file
    if (!gotLine) {
      stop("end of file reached before all expected genotypes were read");
    }
    if (biallelic[i_geno] && phased) {
      // store the phased genotypes, checking that they are all diploid, and
      // that heterozygotes are phased
      int error = 0;
      long bad_indiv = phased_genotype_bytes(line, &alt_counts[0],
                                             &separator_pos[0], n,
                                             (unsigned char) missing_value,
                                             phased_lut_bytes, tab_delimiter,
                                             error);
      if (bad_indiv >= 0) {
        if (is_gz) {
          gzclose(vcfGzFile);
        } else {
          vcfFile.close();
        }
        // describe the locus (chromosome and position) and the individual
        std::vector<std::string> fields = split_locus_meta(line, "\t");
        std::string indiv = (bad_indiv < sample_names.size()) ?
          std::string(sample_names[bad_indiv]) :
          ("number " + std::to_string(bad_indiv + 1));
        std::string where = "individual " + indiv + " at locus " +
          fields[0] + ":" + fields[1];
        if (error == 3) {
          stop("the line of locus " + fields[0] + ":" + fields[1] +
               " has fewer genotypes than individuals");
        } else if (error == 2) {
          stop("the genotype of " + where + " is an unphased heterozygote " +
               "('0/1' or '1/0'): phased genotypes need '|' between the " +
               "alleles of heterozygotes (use gen_tibble() to read " +
               "unphased genotypes)");
        } else {
          stop("the genotype of " + where + " is not diploid: phased " +
               "genotypes can only be read for diploid individuals");
        }
      }
      ptr = std::copy(&alt_counts[0], &alt_counts[0] + n, ptr);
    } else if (biallelic[i_geno]){
      // count alternate alleles, checking that no genotype has more alleles
      // than max_ploidy (in case we estimated it incorrectly in the first line)
      if (!count_alt_alleles(line, &alt_counts[0], &separator_pos[0], n,
                             missing_value, tab_delimiter)) {
        if (is_gz) {
          gzclose(vcfGzFile);
        } else {
          vcfFile.close();
        }
        stop("a genotype has more than max_ploidy alleles. We estimate max_ploidy from the first variant in the vcf file, make sure that variant is representative of ploidy (e.g. it is not on a sex chromosome)");
      }
      ptr = std::copy(&alt_counts[0], &alt_counts[0] + n, ptr);
    }

  }

  // close the file
  if (is_gz) {
    gzclose(vcfGzFile);
  } else {
    vcfFile.close();
  }
  // delete the array
//  delete[] alt_counts;
//  delete[] separator_pos;
  // return true if the file was read successfully
  return true;

}
