#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <climits>
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

// Strings stored one after the other in a single buffer, with the offset of
// the start of each (and of the end of the last). With tens of millions of
// loci names, this takes much less memory than a std::vector<std::string>
// (32 bytes per string, plus a separate allocation for strings longer than
// 15 characters, and copies of all of them whenever the vector grows)
struct StringBuffer {
  std::string chars;
  std::vector<size_t> starts;
  StringBuffer() { starts.push_back(0); }
  void reserve(size_t n_strings, size_t n_chars) {
    starts.reserve(n_strings + 1);
    chars.reserve(n_chars);
  }
  void push_back(const std::string &x) {
    chars.append(x);
    starts.push_back(chars.size());
  }
  void push_back(const char* x, size_t len) {
    chars.append(x, len);
    starts.push_back(chars.size());
  }
  size_t size() const { return starts.size() - 1; }
  // convert to an R character vector (freeing the buffer)
  CharacterVector to_r() {
    CharacterVector res(size());
    for (size_t i = 0; i < size(); i++) {
      SET_STRING_ELT(res, i, Rf_mkCharLen(chars.data() + starts[i],
                                          starts[i + 1] - starts[i]));
    }
    std::string().swap(chars);
    std::vector<size_t>().swap(starts);
    return res;
  }
};

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

// Function to find the metadata of a locus (its first 5 fields: CHROM, POS,
// ID, REF and ALT) in a line, without copying them: it stores the start and
// the length of each field, and returns false if the line has fewer than 5
// fields. As when splitting the line on tabs, consecutive tabs are a single
// separator.
inline bool find_locus_meta(const std::string &line,
                            const char** field_start,
                            size_t* field_len) {
  const char* p = line.c_str();
  const char* line_end = p + line.size();
  // skip leading tabs
  while ((p < line_end) && (*p == '\t')) {
    p++;
  }
  for (int i = 0; i < 5; i++) {
    if (p >= line_end) {
      return false;
    }
    const char* field_end = static_cast<const char*>(
      std::memchr(p, '\t', line_end - p));
    if (field_end == NULL) {
      field_end = line_end;
    }
    field_start[i] = p;
    field_len[i] = field_end - p;
    // move to the start of the next field
    p = field_end;
    while ((p < line_end) && (*p == '\t')) {
      p++;
    }
  }
  return true;
}

// Function to read a position (a 32 bit integer, as in std::stoi), stopping
// with an error if it is not a number
inline int32_t parse_position(const char* p) {
  char* end;
  errno = 0;
  long value = std::strtol(p, &end, 10);
  if ((end == p) || (errno == ERANGE) || (value > INT_MAX) ||
      (value < INT_MIN)) {
    stop("the position of a variant is not a valid integer");
  }
  return static_cast<int32_t>(value);
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
  std::string chromosome_key; // reused to look up a chromosome name
  StringBuffer marker_id;
  if (read_names) {
    marker_id.reserve(10e6, 10e7);
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

    // Find chromosome, position, ID, REF and ALT (the first 5 fields) in the
    // line, without copying them
    const char* field_start[5];
    size_t field_len[5];
    if (!find_locus_meta(line, field_start, field_len)) {
      stop("a variant line of the VCF has fewer than 5 fields");
    }

    // Check if the marker is biallelic SNP (REF and ALT allele are single character)
    // if yes, get the metadata for this locus
    if ((field_len[3] == 1) && (field_len[4] == 1)) {
      biallelic.push_back(true);
      // find the index of the chromosome (loci are usually sorted by
      // chromosome, so we only need a lookup when the chromosome changes)
      if ((this_chromosome < 0) ||
          (chromosome_names[this_chromosome].compare(
              0, std::string::npos, field_start[0], field_len[0]) != 0)) {
        chromosome_key.assign(field_start[0], field_len[0]);
        std::unordered_map<std::string, int>::iterator found =
          chromosome_lookup.find(chromosome_key);
        if (found == chromosome_lookup.end()) {
          this_chromosome = chromosome_names.size();
          chromosome_lookup[chromosome_key] = this_chromosome;
          chromosome_names.push_back(chromosome_key);
        } else {
          this_chromosome = found->second;
        }
      }
      chromosome_index.push_back(this_chromosome + 1);
      physical_pos.push_back(parse_position(field_start[1]));
      if (read_names) {
        marker_id.push_back(field_start[2], field_len[2]);
      }
      allele1.push_back(field_start[4][0]);
      allele2.push_back(field_start[3][0]);
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
    loci_table["marker.ID"] = marker_id.to_r();
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

inline bool count_alt_alleles(const std::string &line,
                              unsigned char* alt_counts,
                              size_t* separator_pos,
                              size_t n,
                              const size_t missingValue,
                              const std::string &separator
) {
  separator_pos[0] = line.find(separator, 0);
  // find the position of all separators
  for (size_t i = 1; i < 8+n; i++) {
    separator_pos[i] = line.find(separator, separator_pos[i-1]+1);
  }
  separator_pos[8+n] = line.length();

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

/**
 * @brief Reads genotype data from a VCF file and stores alternate allele counts in a file-backed matrix.
 *
 * For each biallelic locus in the VCF file, counts the number of alternate alleles per individual and writes these counts to the provided file-backed matrix (FBM). Supports both plain text and gzipped VCF files. Skips the specified number of header lines and processes only loci marked as biallelic. Verifies that no genotype exceeds the maximum ploidy inferred from the first variant line.
 *
 * @param filename Path to the VCF file (plain or gzipped).
 * @param biallelic Logical vector indicating which loci are biallelic and should be processed.
 * @param missing_value Value to assign for missing genotypes.
 * @param n_header_lines Number of header lines to skip in the VCF file.
 * @return true if all genotypes are successfully read and stored.
 *
 * @throws Rcpp::exception if the file cannot be opened, if the end of file is reached prematurely, or if a genotype exceeds the maximum allowed ploidy.
 */

// [[Rcpp::export]]
bool vcf_genotypes_to_fbm(std::string filename,
                          Environment BM,
                          LogicalVector& biallelic,
                          const int missing_value,
                          const int n_header_lines) {
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
    if (biallelic[i_geno]){
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
