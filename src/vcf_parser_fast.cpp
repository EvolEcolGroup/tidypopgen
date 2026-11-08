#include <Rcpp.h>
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
// ID, REF and ALT) in a line of length len, without copying them: it stores
// the start and the length of each field, and returns false if the line has
// fewer than 5 fields. As when splitting the line on tabs, consecutive tabs
// are a single separator.
inline bool find_locus_meta(const char* line, size_t len,
                            const char** field_start,
                            size_t* field_len) {
  const char* p = line;
  const char* line_end = line + len;
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
 * count the number of alternate alleles of each genotype of a line
 *
 * The genotypes start at the 10th field of the line (of length len). Each
 * allele is a single character ('0', '1' or '.'), and alleles are delimited
 * by '/' or '|' (unphased or phased), so we jump from one allele to the next,
 * until we reach the end of the genotype section (':' or the end of the
 * field). A genotype with any missing allele is missing.
 *
 * @param line the line (a C string)
 * @param len the length of the line
 * @param alt_counts the array to store the counts of the n individuals
 * @param n the number of individuals
 * @param missing_value the value stored for missing genotypes
 * @param max_alleles the largest number of alleles of a genotype so far
 *   (updated)
 * @return false if the line has fewer than n genotypes
 */
inline bool count_alt_alleles(const char* line, size_t len,
                              unsigned char* alt_counts, size_t n,
                              unsigned char missing_value,
                              size_t &max_alleles) {
  const char* line_end = line + len;
  // the start of the 10th field
  const char* p = skip_fields(line, 9);
  for (size_t i = 0; i < n; i++) {
    if ((p == NULL) || (p > line_end)) {
      return false;
    }
    const char* field_end = static_cast<const char*>(
      std::memchr(p, '\t', line_end - p));
    if (field_end == NULL) {
      field_end = line_end;
    }
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
    if (n_alleles > max_alleles) {
      max_alleles = n_alleles;
    }
    alt_counts[i] = is_missing ? missing_value : n_alt;
    // move to the next field
    p = field_end + 1;
  }
  return true;
}

// A reader of the lines of a VCF (plain or gzipped), a chunk at a time: each
// chunk holds whole lines, each ended by '\0' (written in place of the
// newline, and of a carriage return before it), so that they can be parsed
// as C strings. The lines of a chunk are valid until the next chunk is read.
class VcfLineReader {
public:
  explicit VcfLineReader(const std::string &filename)
    : fp(NULL), gz(NULL), buffer(1 << 23), carry_start(0), n_carry(0),
      at_eof(false) {
    if (is_gz_file(filename)) {
      gz = gzopen(filename.c_str(), "rb");
      if (gz == NULL) {
        stop("Error opening gzipped file");
      }
      gzbuffer(gz, 1 << 17);
    } else {
      fp = std::fopen(filename.c_str(), "rb");
      if (fp == NULL) {
        stop("Error opening file");
      }
    }
  }
  ~VcfLineReader() {
    if (gz != NULL) {
      gzclose(gz);
    }
    if (fp != NULL) {
      std::fclose(fp);
    }
  }

  // read the next chunk of lines, storing the start and length of each line;
  // returns false at the end of the file
  bool next_chunk(std::vector<char*> &line_start,
                  std::vector<size_t> &line_len) {
    line_start.clear();
    line_len.clear();
    // move the incomplete line left over from the previous chunk to the
    // start of the buffer
    if (n_carry > 0) {
      std::memmove(&buffer[0], &buffer[carry_start], n_carry);
    }
    size_t n_data = n_carry;
    n_carry = 0;
    // read until the buffer holds at least one whole line (growing the
    // buffer if a line is longer than it) or the end of the file
    size_t last_newline = 0;
    bool found_newline = false;
    while (true) {
      if (!at_eof && (n_data < buffer.size())) {
        size_t n_read = read(&buffer[n_data], buffer.size() - n_data);
        if (n_read == 0) {
          at_eof = true;
        }
        n_data += n_read;
      }
      // the last newline in the buffer
      for (size_t i = n_data; i > 0; i--) {
        if (buffer[i - 1] == '\n') {
          last_newline = i - 1;
          found_newline = true;
          break;
        }
      }
      if (found_newline || at_eof) {
        break;
      }
      buffer.resize(buffer.size() * 2);
    }
    if (n_data == 0) {
      return false;
    }
    // the end of the lines of this chunk: after the last newline, or (at the
    // end of the file) the end of the data, after which we add a '\0'
    size_t chunk_end;
    if (at_eof && (!found_newline || (last_newline + 1 < n_data))) {
      if (n_data == buffer.size()) {
        buffer.resize(buffer.size() + 1);
      }
      buffer[n_data] = '\n';
      chunk_end = n_data + 1;
    } else {
      chunk_end = last_newline + 1;
      // keep the incomplete line for the next chunk
      carry_start = chunk_end;
      n_carry = n_data - chunk_end;
    }
    // split the lines, ending each with '\0'
    char* p = &buffer[0];
    char* end = &buffer[0] + chunk_end;
    while (p < end) {
      char* newline = static_cast<char*>(std::memchr(p, '\n', end - p));
      *newline = '\0';
      size_t len = newline - p;
      // remove the carriage return of CRLF line endings
      if ((len > 0) && (p[len - 1] == '\r')) {
        len--;
        p[len] = '\0';
      }
      line_start.push_back(p);
      line_len.push_back(len);
      p = newline + 1;
    }
    return true;
  }

private:
  std::FILE* fp;
  gzFile gz;
  std::vector<char> buffer;
  size_t carry_start;
  size_t n_carry;
  bool at_eof;

  // read up to n bytes into dest, returning the number of bytes read (0 at
  // the end of the file)
  size_t read(char* dest, size_t n) {
    if (gz != NULL) {
      int n_read = gzread(gz, dest, n);
      if (n_read < 0) {
        int err;
        stop(std::string("error reading the gzipped file: ") +
          gzerror(gz, &err));
      }
      return n_read;
    }
    size_t n_read = std::fread(dest, 1, n, fp);
    if ((n_read < n) && std::ferror(fp)) {
      stop("error reading the file");
    }
    return n_read;
  }
};

// A writer that appends bytes to a new file, and removes the file if it is
// destroyed before finish() (e.g. when an error stops the parsing)
class BkWriter {
public:
  explicit BkWriter(const std::string &filename)
    : filename(filename), finished(false) {
    fp = std::fopen(filename.c_str(), "wb");
    if (fp == NULL) {
      stop("could not create the backing file " + filename);
    }
  }
  ~BkWriter() {
    if (!finished) {
      if (fp != NULL) {
        std::fclose(fp);
      }
      std::remove(filename.c_str());
    }
  }
  void write(const std::vector<unsigned char> &x) {
    if ((x.size() > 0) &&
        (std::fwrite(&x[0], 1, x.size(), fp) != x.size())) {
      stop("error writing the backing file " + filename);
    }
  }
  void finish() {
    if (std::fclose(fp) != 0) {
      fp = NULL;
      stop("error writing the backing file " + filename);
    }
    fp = NULL;
    finished = true;
  }

private:
  std::string filename;
  std::FILE* fp;
  bool finished;
};

/******************************************************************************
 * Read a VCF into a loci table and a backing file of genotypes
 *
 * This function reads the VCF once: it creates a loci table with the
 * information on the biallelic snps, finds the sample names and the ploidy of
 * each individual, and writes the number of alternate alleles of each
 * individual at each biallelic snp to the backing file of an FBM (one byte
 * per genotype, the genotypes of a locus one after the other, with max_ploidy
 * + 1 for missing genotypes).
 *
 * The ploidy of each individual is the number of alleles of its genotype in
 * the first variant; for individuals whose genotype is missing (".") there,
 * we use their first genotype that is not "." in the next 1000 variants (or
 * else the largest ploidy of the other individuals, or 2). Until those
 * ploidies are known, the genotypes are kept in memory, with 255 as a
 * provisional value for missing genotypes.
 *
 * To limit memory use with many loci, chromosomes are returned as an integer
 * index (1-based) into a vector of chromosome names (in order of appearance),
 * alleles (a single character for biallelic snps) are stored as characters,
 * and locus names are only read if they are needed.
 *
 * @param filename The name of the VCF file to parse
 * @param bk_file The name of the backing file to create
 * @param read_names Whether to read the names (IDs) of the loci
 */
// [[Rcpp::export]]
List vcf_to_bk(std::string filename, std::string bk_file,
               bool read_names = true) {
  VcfLineReader reader(filename);
  BkWriter writer(bk_file);

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
  std::string last_header_line; // the last line of the header has the samples
  bool first_variant = true;

  // individuals whose genotypes have all been missing (".") so far, and the
  // number of variants in which we still look for their ploidy (we only look
  // in the first few variants, so that the search does not slow down
  // parsing)
  std::vector<size_t> unknown_ploidy;
  int ploidy_lines_left = 1000;
  // whether the ploidies are known (and so the value for missing genotypes)
  bool ploidy_known = false;
  const unsigned char provisional_missing = 255;
  unsigned char missing_value = provisional_missing;
  size_t max_ploidy = 0;
  // largest number of alleles of a genotype, which can not exceed max_ploidy
  size_t max_alleles = 0;

  // genotypes not yet written to the backing file
  std::vector<unsigned char> genotypes;
  const size_t genotypes_flush_size = 1 << 24;
  size_t n_individuals = 0;

  // once the ploidies are known: give the remaining unknown ploidies the
  // largest ploidy of the other individuals (or 2 if no ploidy could be
  // estimated), set the value for missing genotypes, check the genotypes so
  // far, and replace the provisional value for missing genotypes
  auto set_ploidy = [&]() {
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
      unknown_ploidy.clear();
    }
    for (size_t i = 0; i < ploidy.size(); ++i) {
      max_ploidy = std::max(max_ploidy, (size_t) ploidy[i]);
    }
    missing_value = max_ploidy + 1;
    std::replace(genotypes.begin(), genotypes.end(), provisional_missing,
                 missing_value);
    ploidy_known = true;
  };
  // a genotype can not have more alleles than max_ploidy (in case we
  // estimated it incorrectly in the first line)
  auto check_alleles = [&]() {
    if (ploidy_known && (max_alleles > max_ploidy)) {
      stop("a genotype has more than max_ploidy alleles. We estimate max_ploidy from the first variant in the vcf file, make sure that variant is representative of ploidy (e.g. it is not on a sex chromosome)");
    }
  };

  std::vector<char*> line_start;
  std::vector<size_t> line_len;
  while (reader.next_chunk(line_start, line_len)) {
    for (size_t l = 0; l < line_start.size(); l++) {
      char* line = line_start[l];
      size_t len = line_len[l];
      // skip empty lines
      if (len == 0) {
        continue;
      }
      // skip header lines (saving the last one, with the sample names)
      if (line[0] == '#') {
        if (first_variant) {
          last_header_line.assign(line, len);
        }
        continue;
      }

      // Find chromosome, position, ID, REF and ALT (the first 5 fields) in
      // the line, without copying them
      const char* field_start[5];
      size_t field_len[5];
      if (!find_locus_meta(line, len, field_start, field_len)) {
        stop("a variant line of the VCF has fewer than 5 fields");
      }

      if (first_variant) {
        // get the sample names, from the last line of the header
        std::vector<std::string> fields = split_line(last_header_line, "\t");
        for (unsigned int i = 9; i < fields.size(); ++i) {
          sample_names.push_back(fields[i]);
        }
        // get the ploidy of each individual from its genotype (the genotypes
        // start at the 10th field; a missing genotype "." gives no
        // information on ploidy, so for those individuals we use their first
        // genotype that is not ".")
        const char* p = skip_fields(line, 9);
        while (p != NULL) {
          int this_ploidy = genotype_ploidy(p);
          if (this_ploidy == 0) {
            unknown_ploidy.push_back(ploidy.size());
          }
          ploidy.push_back(this_ploidy);
          p = skip_fields(p, 1);
        }
        n_individuals = ploidy.size();
        first_variant = false;
        if (unknown_ploidy.empty()) {
          set_ploidy();
        }
      } else if (!ploidy_known) {
        // look for the ploidy of the individuals with unknown ploidy
        ploidy_lines_left--;
        const char* p = line;
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
        if (unknown_ploidy.empty() || (ploidy_lines_left == 0)) {
          set_ploidy();
          check_alleles();
        }
      }

      // Check if the marker is biallelic SNP (REF and ALT allele are single
      // character); if yes, get the metadata and genotypes of this locus
      if ((field_len[3] == 1) && (field_len[4] == 1)) {
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
        // count the alternate alleles of each individual
        size_t offset = genotypes.size();
        genotypes.resize(offset + n_individuals);
        if (!count_alt_alleles(line, len, &genotypes[offset], n_individuals,
                               missing_value, max_alleles)) {
          stop("a variant line of the VCF has fewer genotypes than individuals");
        }
        check_alleles();
      }
    }
    // write the genotypes to the backing file (once the value for missing
    // genotypes is known)
    if (ploidy_known && (genotypes.size() >= genotypes_flush_size)) {
      writer.write(genotypes);
      genotypes.clear();
    }
    Rcpp::checkUserInterrupt();
  }
  if (!ploidy_known && !first_variant) {
    set_ploidy();
    check_alleles();
  }
  writer.write(genotypes);
  std::vector<unsigned char>().swap(genotypes);
  writer.finish();

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
                       Named("ploidy") = ploidy));
}
