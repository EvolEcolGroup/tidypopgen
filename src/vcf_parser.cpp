#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <zlib.h>

using namespace Rcpp;

// Function to split a string by multiple delimiters and return a vector of tokens
std::vector<std::string> split(const std::string &s, const std::string &delimiters) {
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

// Function to count the number of alternate alleles from genotype information
// missingValue is max_ploidy +1
int countAlternateAlleles(const std::string &genotype, const int missingValue) {
  int count = 0;
  if (genotype[0] == '.'){
    return missingValue;
  }
  // find first : in the string genotype
  int colon_pos = genotype.find_first_of(":");
  //Rcout<<colon_pos<<std::endl;
  if (colon_pos == -1){
    colon_pos = genotype.size();
  }

  //Rcout<<genotype<<std::endl;
  for (int pos_i = 0; pos_i < genotype.size(); pos_i++){
    if (genotype[pos_i] == '1'){
      count++;
    }
    pos_i++;
    //
    if (pos_i< genotype.size()){ // don't read from string if we got beyond the limit
      if (genotype[pos_i] == ':'){
        break;
      }
    }

  }
  return count;
}


// Function to read a line from a gzipped file
bool getline_gz(gzFile file, std::string &line) {
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

// [[Rcpp::export]]
List extractAltAlleleCountsFromVCF(std::string filename,
                                   IntegerMatrix& allele_counts, IntegerVector& ploidy,
                                   int numIndividuals, int missingValue,
                                   int maxLoci = 1000, int skipLoci = 100, bool diploid = false) {

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

  // metadata
  std::vector<std::string>chromosome;
  chromosome.reserve(maxLoci);
  std::vector<std::string> marker_id;
  marker_id.reserve(maxLoci);
  std::vector<std::string> physical_pos;
  physical_pos.reserve(maxLoci);
  std::vector<std::string> allele1;
  allele1.reserve(maxLoci);
  std::vector<std::string> allele2;
  allele2.reserve(maxLoci);

  std::string line;
  int lociCount = 0; // number of biallelic loci
  int lineCount = 0; // number of valid lines
  int skippedLoci = 0;
  bool end_of_file = false;

  // Read the file line by line
  while (lineCount < maxLoci) {
    bool gotLine;
    if (is_gz) {
      gotLine = getline_gz(vcfGzFile, line);
    } else {
      gotLine = (bool)std::getline(vcfFile, line);
    }
    if (!gotLine) {
      end_of_file = true; // we got to the end of the file
      break;
    }

    // Skip header lines
    if (line[0] == '#') {
      continue;
    }

    // Skip the specified number of loci
    if (skippedLoci < skipLoci) {
      skippedLoci++;
      continue;
    }

    // Split the VCF line into fields
    std::vector<std::string> fields = split(line, "\t");

    // Check if the marker is biallelic SNP (REF and ALT allele are single character)
    if ((fields[3].length() ==1) && (fields[4].length() ==1)) {
      chromosome.push_back(fields[0]);
      physical_pos.push_back(fields[1]);
      marker_id.push_back(fields[2]);
      allele1.push_back(fields[4]);
      allele2.push_back(fields[3]);

      // Collect alternate allele counts for all individuals
      if (fields[4][0] != '.'){
        for (unsigned int i = 9; i < fields.size(); ++i) {
          int count = countAlternateAlleles(fields[i], missingValue);
          allele_counts(i - 9, lociCount) = count;
        }
      } else { // if this is an invariant site
        for (unsigned int i = 9; i < fields.size(); ++i) {
          allele_counts(i - 9, lociCount) = 0;
        }
      }


      lociCount++;
    }
    lineCount++;
  }
  // test if the next line would have been the end of file
  bool gotLine;
  if (is_gz) {
    gotLine = getline_gz(vcfGzFile, line);
  } else {
    gotLine = (bool)std::getline(vcfFile, line);
  }
  if (!gotLine) {
    end_of_file = true; // we got to the end of the file
  }

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

  // Return a list with useful information
  return List::create(
    Named("loci_table") = loci_table,
    Named("num_loci") = lociCount,
    Named("end_of_file") = end_of_file
  );
}


// get ploid of each individual from the first line of data
// [[Rcpp::export]]
List get_ploidy_from_VCF(std::string filename) {
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

  std::string line;
  // when we find the genotypes, we need to look at the old line for the names of samples
  std::string line_old;
  int numIndividuals = 0;
  std::vector<int> ploidy;
  std::vector<std::string> sample_names;

  // Read the file line by line
  while (true) {
    bool gotLine;
    if (is_gz) {
      gotLine = getline_gz(vcfGzFile, line);
    } else {
      gotLine = (bool)std::getline(vcfFile, line);
    }
    if (!gotLine) {
      break;
    }

    // Skip header lines
    if (line[0] == '#') {
      line_old = line;
      continue;
    }

    // Split the VCF line into fields
    std::vector<std::string> fields = split(line, "\t");
    numIndividuals = fields.size() - 9; // 9 fixed fields before individual data
    ploidy.resize(numIndividuals, 0);
    for (unsigned int i = 9; i < fields.size(); ++i) {
      ploidy[i - 9] = split(fields[i], "/|").size();
      //ploidy[i - 9] = 1+ (fields[i].length()-1)/2;
    }
    // now process names
    fields = split(line_old, "\t");
    sample_names.resize(numIndividuals,"");
    for (unsigned int i = 9; i < fields.size(); ++i) {
      sample_names[i - 9] = fields[i];
    }
    break;
  }

  if (is_gz) {
    gzclose(vcfGzFile);
  } else {
    vcfFile.close();
  }

  // Return a list containing the matrix and the ploidy vector
  return (List::create(Named("sample_names") = sample_names,
                       Named("ploidy") = ploidy));
}

