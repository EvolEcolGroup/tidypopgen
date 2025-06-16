#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
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
 * @param filename The name of the VCF file to parse
 */
// [[Rcpp::export]]
List vcf_loci_table(std::string filename) {
  // open the file
  bool is_gz = (filename.substr(filename.size() - 3) == ".gz");
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
  // loci metadata (preallocate a million loci, as a start)
  std::vector<std::string>chromosome;
  chromosome.reserve(10e6);
  std::vector<std::string> marker_id;
  marker_id.reserve(10e6);
  std::vector<int32_t> physical_pos;
  physical_pos.reserve(10e6);
  std::vector<std::string> allele1;
  allele1.reserve(10e6);
  std::vector<std::string> allele2;
  allele2.reserve(10e6);
  // individual metadata
  std::vector<int> ploidy;
  std::vector<std::string> sample_names;
  
  
  std::string line;
  std::string line_old; // used to get the samples from the last line of the header
  int n_header_lines =0;
  
  bool get_sample_names = true;
  bool get_ploidy = true;
  
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
      chromosome.push_back(fields[0]);
      physical_pos.push_back(std::stoi(fields[1]));
      marker_id.push_back(fields[2]);
      allele1.push_back(fields[4]);
      allele2.push_back(fields[3]);
    } else {
      biallelic.push_back(false);
    }
    
    
    // Check if we need to get ploidy (this is only the case for the first line
    // of real data)
    if (get_ploidy) {
      std::vector<std::string> fields = split_line(line, "\t");
      for (unsigned int i = 9; i < fields.size(); ++i) {
        ploidy.push_back(split_line(fields[i], "/|").size());
      }
      get_ploidy = false;
    }
    
  }
  
  // close the file
  if (is_gz) {
    gzclose(vcfGzFile);
  } else {
    vcfFile.close();
  }
  
  // craete loci table
  DataFrame loci_table = DataFrame::create( Named("chromosome") = chromosome ,
                                            _["marker.ID"] = marker_id,
                                            _["physical.pos"] = physical_pos,
                                            _["allele1"] = allele1,
                                            _["allele2"] = allele2);
  
  
  // Return a list containing the matrix and the ploidy vector
  return (List::create(Named("loci_tbl") = loci_table,
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
 */

inline void count_alt_alleles(const std::string &line,
                              unsigned char* alt_counts,
                              size_t* separator_pos,
                              size_t n,
                              const size_t missingValue,
                              const std::string &separator
) {
  // reset counts to 0
  std::fill(alt_counts, alt_counts+n, 0);
  // for (uint i = 0; i < n; i++) {
  //   alt_counts[i] = 0;
  // }
  
  separator_pos[0] = line.find(separator, 0);
  // find the position of all separators
  for (size_t i = 1; i < 8+n; i++) {
    separator_pos[i] = line.find(separator, separator_pos[i-1]+1);
  }
  separator_pos[8+n] = line.length();
  
  // now move through each of the genotypes and process them
  // we ignore the first 9 fields (so start from position 9)
  for (size_t i=0; i<n; i++){
    size_t pos = separator_pos[i+8];
    // the genotypes start at pos +1
    if (line[pos+1] == '.') {
      alt_counts[i] = missingValue;
    } else {
      for (size_t pos_i = pos+1; pos_i < separator_pos[i+9]+1; pos_i++){
        if (line[pos_i] == '1'){
          alt_counts[i]++;
        } else if (line[pos_i] == ':') //if at the end of the genotype section
          break;
      }
    }
  }
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
                          IntegerVector& biallelic,
                          const int missing_value,
                          const int n_header_lines) {
  // open the file
  bool is_gz = (filename.substr(filename.size() - 3) == ".gz");
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
    if (is_gz) {
      gotLine = gz_getline(vcfGzFile, line);
    } else {
      gotLine = (bool)std::getline(vcfFile, line);
    }
    // Check if we reached the end of the file
    if (!gotLine) {
      stop("end of file reached before all expected genotypes were read");
    }
    if (biallelic[i_geno]){
      count_alt_alleles(line, &alt_counts[0], &separator_pos[0], n,missing_value, tab_delimiter);
      ptr = std::copy(&alt_counts[0], &alt_counts[0] + n, ptr);
    }
    
  }
  // check ploidy for the last line, in case we estimated it incorrectly in the first
  // line
  std::vector<std::string> fields = split_line(line, "\t");
  for (unsigned int i = 9; i < fields.size(); ++i) {
    if (split_line(fields[i], "/|").size() > size_t((missing_value-1))) {
      stop("a genotype has more than max_ploidy alleles. We estimate max_ploidy from the first variant in the vcf file, make sure that variant is representative of ploidy (e.g. it is not on a sex chromosome)");
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
