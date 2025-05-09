#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <zlib.h>

using namespace Rcpp;

// Function to read a line from a gzipped file
bool gz_getline(gzFile file, std::string &line) {
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
std::vector<std::string> split_line(const std::string &s, const std::string &delimiters) {
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
std::vector<std::string> split_locus_meta(const std::string &s, const std::string &delimiters) {
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
      // next line
      continue;
    }
    
    // Extract chromosome, position, and 
    // Split the VCF line into fields (only for the first few fields)
    std::vector<std::string> fields = split_locus_meta(line, "\t");
    
    // Check if the marker is biallelic SNP (REF and ALT allele are single character)
    // if yes, get the metadata for this locus
    if ((fields[3].length() ==1) && (fields[4].length() ==1)) {
      chromosome.push_back(fields[0]);
      physical_pos.push_back(std::stoi(fields[1]));
      marker_id.push_back(fields[2]);
      allele1.push_back(fields[4]);
      allele2.push_back(fields[3]);
    }
    
    
    // Check if we need to get ploidy (this is only the case for the first line
    // of real data)
    if (get_ploidy) {
      std::vector<std::string> fields = split_line(line, "\t");
      int numIndividuals = fields.size() - 9; // 9 fixed fields before individual data
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
  return (List::create(Named("sample_names") = sample_names,
                       Named("ploidy") = ploidy,
                       Named("loci_tbl") = loci_table));
}

