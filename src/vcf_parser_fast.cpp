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
#include <cstdint>
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

// Function to check whether the first 18 bytes of a file are the header of a
// BGZF block (the blocked gzip format of bgzip and htslib, as checked by
// htslib): a gzip member with a 6-byte extra field "BC" holding the size of
// the block
inline bool is_bgzf_header(const unsigned char* h) {
  return (h[0] == 31) && (h[1] == 139) && (h[2] == 8) && ((h[3] & 4) != 0) &&
    (h[10] == 6) && (h[11] == 0) && (h[12] == 'B') && (h[13] == 'C') &&
    (h[14] == 2) && (h[15] == 0);
}

// Function to read a little-endian unsigned integer of n bytes
inline size_t read_le(const unsigned char* p, int n) {
  size_t x = 0;
  for (int i = n - 1; i >= 0; i--) {
    x = (x << 8) | p[i];
  }
  return x;
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

// The names of the chromosomes, in order of appearance, with a hash table to
// find the index of a name. The names are kept in a StringBuffer, and the
// hash table (open addressing, with linear probing) only holds their indices,
// so that nothing is allocated for each chromosome: with millions of
// chromosomes (e.g. RAD loci assembled de novo, one per locus), allocating
// (and then freeing) a std::string and a node of a std::unordered_map for
// each took much of the time (and memory) of reading the VCF
class ChromosomeTable {
public:
  StringBuffer names;
  ChromosomeTable() : slots(1 << 10, -1) {}
  // whether the name of chromosome k is the string name (of length len)
  bool is(int k, const char* name, size_t len) const {
    return (names.starts[k + 1] - names.starts[k] == len) &&
      (std::memcmp(names.chars.data() + names.starts[k], name, len) == 0);
  }
  // the index of the chromosome called name, adding it if it is new
  int find_or_add(const char* name, size_t len) {
    uint32_t h = hash(name, len);
    size_t mask = slots.size() - 1;
    size_t i = h & mask;
    while (slots[i] >= 0) {
      int k = slots[i];
      if ((hashes[k] == h) && is(k, name, len)) {
        return k;
      }
      i = (i + 1) & mask;
    }
    int k = names.size();
    names.push_back(name, len);
    hashes.push_back(h);
    slots[i] = k;
    // keep the table less than 70% full (so that probes stay short)
    if (10 * hashes.size() > 7 * slots.size()) {
      grow();
    }
    return k;
  }

private:
  std::vector<int> slots; // the chromosome in each slot (-1 if empty)
  std::vector<uint32_t> hashes; // the hash of the name of each chromosome
  // FNV-1a hash of a string, with the final mix of MurmurHash3 (so that all
  // its bits, and so the slot, depend on all the characters)
  static uint32_t hash(const char* s, size_t len) {
    uint64_t h = 14695981039346656037ULL;
    for (size_t i = 0; i < len; i++) {
      h = (h ^ (unsigned char) s[i]) * 1099511628211ULL;
    }
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    return (uint32_t) h;
  }
  // double the size of the hash table, and put back the chromosomes
  void grow() {
    slots.assign(slots.size() * 2, -1);
    size_t mask = slots.size() - 1;
    for (size_t k = 0; k < hashes.size(); k++) {
      size_t i = hashes[k] & mask;
      while (slots[i] >= 0) {
        i = (i + 1) & mask;
      }
      slots[i] = k;
    }
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

// Function to check whether a variant is a biallelic SNP (its REF and ALT
// alleles are a single character each), given the lengths of its first 5
// fields (as found by find_locus_meta()). Only these variants are read; the
// two passes of gt_merge_vcfs() use the same test, so that they number the
// SNPs of a file in the same way
inline bool is_biallelic_snp(const size_t* field_len) {
  return (field_len[3] == 1) && (field_len[4] == 1);
}

// Function to read a position (a 32 bit integer, as in std::stoi) from its
// field (of length len), stopping with an error if the field is not a number
// (e.g. "12junk", which would otherwise be read as 12)
inline int32_t parse_position(const char* p, size_t len) {
  char* end;
  errno = 0;
  long value = std::strtol(p, &end, 10);
  if ((end == p) || (end != p + len) || (errno == ERANGE) ||
      (value > INT_MAX) || (value < INT_MIN)) {
    stop("the position of a variant is not a valid integer");
  }
  return static_cast<int32_t>(value);
}

// the outcome of count_alt_alleles()
enum GenotypeStatus {
  GENOTYPES_OK = 0,
  // the line has fewer genotypes than individuals
  TOO_FEW_GENOTYPES = 1,
  // a genotype is not valid for a biallelic SNP (e.g. "2/1", "a/1", "0/12",
  // or an empty field)
  INVALID_GENOTYPE = 2,
  // the FORMAT of the line does not start with GT
  NO_GT_FORMAT = 3
};

/******************************************************************************
 * count the number of alternate alleles of each genotype of a line
 *
 * The genotypes start at the 10th field of the line (of length len), and the
 * genotype (GT) must be the first value of each field (as the VCF
 * specification requires, which we check in the FORMAT field). Each allele
 * of a biallelic SNP is a single character ('0', '1' or '.'), and alleles are
 * delimited by '/' or '|' (unphased or phased), so we jump from one allele to
 * the next, until we reach the end of the genotype section (':' or the end
 * of the field). A genotype with any missing allele is missing.
 *
 * @param line the line (a C string)
 * @param len the length of the line
 * @param alt_counts the array to store the counts of the n individuals
 * @param n the number of individuals
 * @param missing_value the value stored for missing genotypes
 * @param alt_allele the allele that is counted: '1' (the ALT allele), or
 *   '0' to count the REF allele (when the alleles of the SNP are swapped to
 *   match another file, the REF allele becomes the alternate allele)
 * @param max_alleles the largest number of alleles of a genotype so far
 *   (updated)
 * @return a GenotypeStatus (GENOTYPES_OK if the genotypes are valid)
 */
inline int count_alt_alleles(const char* line, size_t len,
                             unsigned char* alt_counts, size_t n,
                             unsigned char missing_value, char alt_allele,
                             size_t &max_alleles) {
  if (n == 0) {
    return GENOTYPES_OK;
  }
  const char* line_end = line + len;
  // the FORMAT field (the 9th) must start with GT
  const char* p = skip_fields(line, 8);
  if ((p == NULL) || (p >= line_end)) {
    return TOO_FEW_GENOTYPES;
  }
  if ((p[0] != 'G') || (p[1] != 'T') ||
      ((p[2] != ':') && (p[2] != '\t') && (p[2] != '\0'))) {
    return NO_GT_FORMAT;
  }
  // the start of the 10th field
  p = skip_fields(p, 1);
  bool invalid = false;
  for (size_t i = 0; i < n; i++) {
    if ((p == NULL) || (p > line_end)) {
      return TOO_FEW_GENOTYPES;
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
      char allele = p[0];
      n_alt += (allele == alt_allele);
      // if any allele is missing, the genotype is missing
      is_missing |= (allele == '.');
      // any other allele (or an empty genotype) is not valid
      invalid |= (allele != '0') & (allele != '1') & (allele != '.');
      n_alleles++;
      if ((p + 1 >= field_end) || ((p[1] != '/') && (p[1] != '|'))) {
        break;
      }
      p += 2;
    }
    // the genotype ends with its field, or with ':' (an allele has a single
    // character, e.g. not "0/12")
    invalid |= (p + 1 < field_end) && (p[1] != ':');
    if (n_alleles > max_alleles) {
      max_alleles = n_alleles;
    }
    alt_counts[i] = is_missing ? missing_value : n_alt;
    // move to the next field
    p = field_end + 1;
  }
  return invalid ? INVALID_GENOTYPE : GENOTYPES_OK;
}

// Function to stop with the error of a GenotypeStatus (other than
// GENOTYPES_OK), for the VCF `what` (e.g. "the VCF" or "the VCF <file>")
inline void stop_genotype_status(int status, const std::string &what) {
  if (status == TOO_FEW_GENOTYPES) {
    stop("a variant line of " + what +
         " has fewer genotypes than individuals");
  }
  if (status == INVALID_GENOTYPE) {
    stop("a genotype of " + what + " is not valid: the alleles of a " +
         "biallelic SNP should be 0, 1 or . (missing), separated by / or |");
  }
  if (status == NO_GT_FORMAT) {
    stop("the FORMAT of a variant of " + what + " does not start with GT " +
         "(the genotype should be the first value of each sample)");
  }
}

// A reader of the lines of a VCF (plain or gzipped), a chunk at a time: each
// chunk holds whole lines, each ended by '\0' (written in place of the
// newline, and of a carriage return before it), so that they can be parsed
// as C strings. The lines of a chunk are valid until the next chunk is read.
// Files compressed with bgzip (BGZF, made of independent blocks of at most
// 64 KB) are read a batch of blocks at a time, and the blocks of a batch are
// decompressed in parallel (with n_threads threads); other gzipped files are
// decompressed sequentially with zlib.
class VcfLineReader {
public:
  VcfLineReader(const std::string &filename, int n_threads)
    : fp(NULL), gz(NULL), is_bgzf(false), n_threads(n_threads),
      buffer(chunk_bytes), carry_start(0), n_carry(0), at_eof(false) {
    // (a thread count that is not positive, e.g. NA, would make OpenMP
    // fail, or try to start a huge number of threads)
    if (n_threads < 1) {
      stop("'n_threads' should be a positive integer");
    }
    if (is_gz_file(filename)) {
      // check whether the file is in BGZF format
      fp = std::fopen(filename.c_str(), "rb");
      if (fp == NULL) {
        stop("Error opening gzipped file");
      }
      unsigned char header[18];
      is_bgzf = (std::fread(header, 1, 18, fp) == 18) &&
        is_bgzf_header(header);
      if (is_bgzf) {
        std::rewind(fp);
      } else {
        std::fclose(fp);
        fp = NULL;
        gz = gzopen(filename.c_str(), "rb");
        if (gz == NULL) {
          stop("Error opening gzipped file");
        }
        gzbuffer(gz, 1 << 17);
      }
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
      if (!at_eof && is_bgzf) {
        // (at least 1 MB, if the buffer is already full of a long line)
        size_t min_bytes = (n_data < chunk_bytes) ? chunk_bytes - n_data :
          (1 << 20);
        n_data += read_bgzf(n_data, min_bytes);
      } else if (!at_eof && (n_data < buffer.size())) {
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
      // (a batch of BGZF blocks grows the buffer as needed)
      if (!is_bgzf) {
        buffer.resize(buffer.size() * 2);
      }
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
  static const size_t chunk_bytes = 1 << 23;
  std::FILE* fp;
  gzFile gz;
  bool is_bgzf;
  int n_threads;
  std::vector<char> buffer;
  size_t carry_start;
  size_t n_carry;
  bool at_eof;
  // BGZF blocks of the current batch (compressed), and the offset of the
  // uncompressed data of each in the batch
  std::vector<std::vector<unsigned char> > blocks;
  std::vector<size_t> block_offset;

  // read the next BGZF block into block (returning false at the end of the
  // file), checking its header
  bool read_bgzf_block(std::vector<unsigned char> &block) {
    unsigned char header[18];
    size_t n_read = std::fread(header, 1, 18, fp);
    if (n_read == 0) {
      if (std::ferror(fp)) {
        stop("error reading the file");
      }
      return false;
    }
    if ((n_read < 18) || !is_bgzf_header(header)) {
      stop("the gzipped file is not a valid BGZF file (as made by bgzip)");
    }
    size_t block_size = read_le(header + 16, 2) + 1;
    if (block_size < 26) {
      stop("the gzipped file is not a valid BGZF file (as made by bgzip)");
    }
    block.resize(block_size);
    std::memcpy(&block[0], header, 18);
    if (std::fread(&block[18], 1, block_size - 18, fp) != block_size - 18) {
      stop("the BGZF file is truncated");
    }
    return true;
  }

  // read BGZF blocks with at least min_bytes of uncompressed data (or until
  // the end of the file), and decompress them in parallel, appending them to
  // the buffer at n_data; returns the number of bytes appended.
  // A possible speed-up, not used for now to avoid a dependency: libdeflate
  // (https://github.com/ebiggers/libdeflate; R packages can link to the
  // static library built by the CRAN package libdeflate, located with
  // system.file() in a configure script) inflates BGZF blocks about 2.2
  // times as fast as zlib, and computes their CRC32 about 3.5 times as fast
  // (measured on a 3.4 GB VCF, in November 2026): reading a bgzipped VCF of
  // 2.3 million variants and 56 samples would take about 3 s less on one
  // core (16.8 s), and about 1 s less on four (12.4 s). libdeflate only
  // decompresses whole buffers of known size, which suits BGZF blocks (their
  // size is in their header, and the uncompressed size in their footer), but
  // not other gzipped files, which would still be read with zlib.
  size_t read_bgzf(size_t n_data, size_t min_bytes) {
    size_t n_blocks = 0;
    size_t total = 0;
    block_offset.assign(1, 0);
    while (total < min_bytes) {
      if (n_blocks == blocks.size()) {
        blocks.resize(n_blocks + 64);
      }
      if (!read_bgzf_block(blocks[n_blocks])) {
        at_eof = true;
        break;
      }
      // the uncompressed size of the block is in its last 4 bytes (at most
      // 64 KB; a larger size, from a corrupted block, would make us allocate
      // up to 4 GB before failing to decompress it)
      const std::vector<unsigned char> &block = blocks[n_blocks];
      size_t block_output_size = read_le(&block[block.size() - 4], 4);
      if (block_output_size > 65536) {
        stop("error decompressing the BGZF file (it may be corrupted)");
      }
      total += block_output_size;
      block_offset.push_back(total);
      n_blocks++;
    }
    if (buffer.size() < n_data + total) {
      buffer.resize(n_data + total);
    }
    // decompress the blocks (raw deflate data between the 18-byte header
    // and the 8-byte footer, CRC32 and size), checking their CRC32
    std::vector<char> block_ok(n_blocks, 1);
    char* out = &buffer[n_data];
#pragma omp parallel num_threads(n_threads)
    {
      z_stream zs;
      std::memset(&zs, 0, sizeof(zs));
      bool zs_ok = (inflateInit2(&zs, -15) == Z_OK);
#pragma omp for schedule(dynamic, 4)
      for (long i = 0; i < (long) n_blocks; i++) {
        std::vector<unsigned char> &block = blocks[i];
        size_t out_size = block_offset[i + 1] - block_offset[i];
        Bytef* this_out = reinterpret_cast<Bytef*>(out + block_offset[i]);
        if (out_size > 0) {
          if (!zs_ok || (inflateReset(&zs) != Z_OK)) {
            block_ok[i] = 0;
            continue;
          }
          zs.next_in = &block[18];
          zs.avail_in = block.size() - 26;
          zs.next_out = this_out;
          zs.avail_out = out_size;
          if ((inflate(&zs, Z_FINISH) != Z_STREAM_END) ||
              (zs.avail_out != 0)) {
            block_ok[i] = 0;
            continue;
          }
        }
        uLong crc = crc32(0L, this_out, out_size);
        if (crc != read_le(&block[block.size() - 8], 4)) {
          block_ok[i] = 0;
        }
      }
      if (zs_ok) {
        inflateEnd(&zs);
      }
    }
    for (size_t i = 0; i < n_blocks; i++) {
      if (!block_ok[i]) {
        stop("error decompressing the BGZF file (it may be corrupted)");
      }
    }
    return total;
  }

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
 * The lines are read a chunk at a time: the metadata of the loci are read in
 * order, and the genotypes of the biallelic snps of a chunk are then counted
 * in parallel (with n_threads threads), as are the blocks of a BGZF file
 * decompressed.
 *
 * With read_genotypes = false, only the loci, the samples and their ploidy
 * are read (no backing file is created): this is the first pass of
 * gt_merge_vcfs() and gt_concat_vcfs(), whose second pass
 * (vcf_merge_genotypes()) reads the genotypes.
 *
 * @param filename The name of the VCF file to parse
 * @param bk_file The name of the backing file to create (ignored if
 *   read_genotypes is false)
 * @param read_names Whether to read the names (IDs) of the loci
 * @param n_threads The number of threads
 * @param read_genotypes Whether to read the genotypes (into bk_file)
 * @return a list with the loci table (`loci_tbl`), `chromosome_names`,
 *   `sample_names`, the `ploidy` of each individual, and `ploidy_guessed`,
 *   whether its ploidy was guessed (its genotypes were all missing in the
 *   first variants)
 */
// [[Rcpp::export]]
List vcf_to_bk(std::string filename, std::string bk_file,
               bool read_names = true, int n_threads = 1,
               bool read_genotypes = true) {
  VcfLineReader reader(filename, n_threads);
  // (the backing file is only created if we read the genotypes)
  std::unique_ptr<BkWriter> writer;
  if (read_genotypes) {
    writer.reset(new BkWriter(bk_file));
  }

  // loci metadata (preallocate 10 million loci, as a start)
  // chromosomes are stored as an index into the table of their names
  std::vector<int> chromosome_index;
  chromosome_index.reserve(10e6);
  ChromosomeTable chromosomes;
  int this_chromosome = -1; // index of the chromosome of the previous locus
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
  // individual metadata (and whether the ploidy of each individual was
  // guessed, because its genotypes were all missing in the first variants)
  std::vector<int> ploidy;
  std::vector<int> ploidy_guessed;
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
    ploidy_guessed.assign(ploidy.size(), 0);
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
          ploidy_guessed[i] = 1;
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
  // the lines of a chunk that are biallelic snps
  std::vector<size_t> snp_lines;
  while (reader.next_chunk(line_start, line_len)) {
    snp_lines.clear();
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
        // each sample in the header needs a genotype (otherwise the error
        // would only come once the FBM is created, from its number of rows)
        if (n_individuals != sample_names.size()) {
          stop("the first variant line of the VCF has " +
               std::to_string(n_individuals) + " genotypes, but the header " +
               "line (#CHROM) has " + std::to_string(sample_names.size()) +
               " samples");
        }
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
      if (is_biallelic_snp(field_len)) {
        // find the index of the chromosome (loci are usually sorted by
        // chromosome, so we only need a lookup when the chromosome changes)
        if ((this_chromosome < 0) ||
            !chromosomes.is(this_chromosome, field_start[0], field_len[0])) {
          this_chromosome = chromosomes.find_or_add(field_start[0],
                                                    field_len[0]);
        }
        chromosome_index.push_back(this_chromosome + 1);
        physical_pos.push_back(parse_position(field_start[1], field_len[1]));
        if (read_names) {
          marker_id.push_back(field_start[2], field_len[2]);
        }
        allele1.push_back(field_start[4][0]);
        allele2.push_back(field_start[3][0]);
        snp_lines.push_back(l);
      }
    }
    if (!read_genotypes) {
      Rcpp::checkUserInterrupt();
      continue;
    }
    // count the alternate alleles of each individual at the biallelic snps
    // of this chunk, in parallel (with the provisional value for missing
    // genotypes if the ploidies are not yet known)
    size_t offset = genotypes.size();
    genotypes.resize(offset + snp_lines.size() * n_individuals);
    unsigned char this_missing = ploidy_known ? missing_value :
      provisional_missing;
    long n_snp_lines = snp_lines.size();
    size_t chunk_max_alleles = 0;
    int status = GENOTYPES_OK;
#pragma omp parallel for num_threads(n_threads) schedule(static) \
    reduction(max:chunk_max_alleles) reduction(max:status)
    for (long k = 0; k < n_snp_lines; k++) {
      size_t l = snp_lines[k];
      int line_status = count_alt_alleles(
        line_start[l], line_len[l], &genotypes[offset + k * n_individuals],
        n_individuals, this_missing, '1', chunk_max_alleles);
      status = std::max(status, line_status);
    }
    if (status != GENOTYPES_OK) {
      stop_genotype_status(status, "the VCF");
    }
    max_alleles = std::max(max_alleles, chunk_max_alleles);
    check_alleles();
    // write the genotypes to the backing file (once the value for missing
    // genotypes is known)
    if (ploidy_known && (genotypes.size() >= genotypes_flush_size)) {
      writer->write(genotypes);
      genotypes.clear();
    }
    Rcpp::checkUserInterrupt();
  }
  if (!ploidy_known && !first_variant) {
    set_ploidy();
    check_alleles();
  }
  if (read_genotypes) {
    writer->write(genotypes);
    std::vector<unsigned char>().swap(genotypes);
    writer->finish();
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

  // Return a list containing the loci table and the ploidy vector (with
  // whether each ploidy was guessed)
  return (List::create(Named("loci_tbl") = loci_table,
                       Named("chromosome_names") = chromosomes.names.to_r(),
                       Named("sample_names") = sample_names,
                       Named("ploidy") = ploidy,
                       Named("ploidy_guessed") = LogicalVector(
                         ploidy_guessed.begin(), ploidy_guessed.end())));
}

/******************************************************************************
 * Write the genotypes of a VCF into the FBM of several merged VCFs
 *
 * This is the second pass of gt_merge_vcfs() and gt_concat_vcfs(): their
 * first pass read the loci of each VCF (with vcf_to_bk(read_genotypes =
 * false)) and matched them across files, giving the column of the merged FBM
 * of each biallelic snp of this file. This function reads the VCF again, and
 * writes the number of alternate alleles of each of its individuals at each
 * of its biallelic snps that is kept straight into its column of the merged
 * FBM, in the rows of this file (row_offset onwards, or the rows given by
 * `row`). The FBM is column-major, so each snp is a contiguous block of
 * n_individuals bytes (unless the rows are given by `row`, in which case the
 * genotypes of a snp are counted in a buffer, then copied to their rows); the
 * snps of a chunk of lines go to different columns, so they are written in
 * parallel (with n_threads threads), as are the blocks of a BGZF file
 * decompressed.
 *
 * The snps are numbered as in the first pass (in the order of the file, with
 * the same test, is_biallelic_snp()), so the file must not change between
 * the two passes (we check that it has the same number of biallelic snps).
 *
 * @param filename The name of the VCF file to read
 * @param BM The merged FBM.code256 (with at least row_offset +
 *   length(ploidy) rows)
 * @param row_offset The (0-based) row of the first individual of this file
 *   (ignored if `row` is not empty)
 * @param row The (0-based) row of each individual of this file in the merged
 *   FBM (e.g. to put the individuals of the file in the order of the FBM), or
 *   an empty vector (rows row_offset onwards)
 * @param col The (0-based) column of each biallelic snp of the file in the
 *   merged FBM, or -1 for snps that are not kept
 * @param swap Whether the alleles of each snp are swapped to match the merged
 *   locus (its REF alleles are counted, rather than its ALT alleles), or an
 *   empty vector (no swapping)
 * @param ploidy The ploidy of each individual of the file (from the first
 *   pass)
 * @param missing_value The value of missing genotypes in the merged FBM
 *   (its largest ploidy + 1)
 * @param n_threads The number of threads
 */
// [[Rcpp::export]]
void vcf_merge_genotypes(std::string filename, Environment BM, int row_offset,
                         const IntegerVector& row, const IntegerVector& col,
                         const LogicalVector& swap,
                         const IntegerVector& ploidy, int missing_value,
                         int n_threads = 1) {
  XPtr<FBM_RW> xp_bm = BM["address_rw"];
  unsigned char* bm = static_cast<unsigned char*>(xp_bm->matrix());
  size_t n_rows = xp_bm->nrow(); // number of rows of the merged FBM
  size_t n_cols = xp_bm->ncol(); // number of columns of the merged FBM
  size_t n_individuals = ploidy.size();
  long n_snps = col.size();
  bool any_swap = swap.size() > 0;
  bool reorder = row.size() > 0;

  // check that the rows and columns of this file are within the FBM
  std::vector<size_t> row_vec;
  if (reorder) {
    if ((size_t) row.size() != n_individuals) {
      stop("'row' should have one value per individual");
    }
    row_vec.resize(n_individuals);
    for (size_t i = 0; i < n_individuals; i++) {
      if ((row[i] == NA_INTEGER) || (row[i] < 0) ||
          ((size_t) row[i] >= n_rows)) {
        stop("'row' is out of the bounds of the merged FBM");
      }
      row_vec[i] = row[i];
    }
  } else if ((row_offset < 0) ||
             ((size_t) row_offset + n_individuals > n_rows)) {
    stop("the merged FBM does not have enough rows for this file");
  }
  for (long k = 0; k < n_snps; k++) {
    if ((col[k] != NA_INTEGER) && (col[k] >= 0) &&
        ((size_t) col[k] >= n_cols)) {
      stop("'col' is out of the bounds of the merged FBM");
    }
  }
  if (any_swap && (swap.size() != n_snps)) {
    stop("'swap' should have one value per biallelic snp");
  }
  if ((missing_value < 1) || (missing_value > 255)) {
    stop("'missing_value' should be between 1 and 255");
  }
  // a genotype can not have more alleles than the largest ploidy of the file
  // (as checked by vcf_to_bk())
  size_t max_ploidy = 0;
  for (size_t i = 0; i < n_individuals; i++) {
    max_ploidy = std::max(max_ploidy, (size_t) ploidy[i]);
  }
  // (copies, so that the parallel loop does not touch R objects)
  std::vector<int> col_vec(col.begin(), col.end());
  std::vector<char> swap_vec(n_snps, 0);
  if (any_swap) {
    for (long k = 0; k < n_snps; k++) {
      swap_vec[k] = (swap[k] == TRUE);
    }
  }
  const unsigned char missing = (unsigned char) missing_value;

  VcfLineReader reader(filename, n_threads);
  std::vector<char*> line_start;
  std::vector<size_t> line_len;
  // the lines of a chunk that are biallelic snps kept in the merge, and the
  // number (in the file) of each of these snps
  std::vector<size_t> snp_lines;
  std::vector<long> snp_number;
  // with reordered rows, the genotypes of the snps of a chunk (one slice of
  // n_individuals per snp, so that the threads write to different slices)
  std::vector<unsigned char> scratch;
  long k = 0; // number of biallelic snps read so far
  while (reader.next_chunk(line_start, line_len)) {
    snp_lines.clear();
    snp_number.clear();
    for (size_t l = 0; l < line_start.size(); l++) {
      // skip empty lines and header lines
      if ((line_len[l] == 0) || (line_start[l][0] == '#')) {
        continue;
      }
      const char* field_start[5];
      size_t field_len[5];
      if (!find_locus_meta(line_start[l], line_len[l], field_start,
                           field_len)) {
        stop("a variant line of the VCF has fewer than 5 fields");
      }
      if (!is_biallelic_snp(field_len)) {
        continue;
      }
      if (k >= n_snps) {
        stop("the VCF " + filename + " has more biallelic SNPs than when " +
             "its loci were read (it changed while being merged)");
      }
      if ((col_vec[k] != NA_INTEGER) && (col_vec[k] >= 0)) {
        snp_lines.push_back(l);
        snp_number.push_back(k);
      }
      k++;
    }
    // count the alternate alleles of each individual at the kept snps of
    // this chunk, in parallel, straight into their columns of the FBM (or,
    // with reordered rows, into the scratch buffer, and then into their rows)
    long n_snp_lines = snp_lines.size();
    if (reorder) {
      scratch.resize((size_t) n_snp_lines * n_individuals);
    }
    size_t chunk_max_alleles = 0;
    int status = GENOTYPES_OK;
#pragma omp parallel for num_threads(n_threads) schedule(static) \
    reduction(max:chunk_max_alleles) reduction(max:status)
    for (long s = 0; s < n_snp_lines; s++) {
      size_t l = snp_lines[s];
      long snp = snp_number[s];
      unsigned char* col_start = bm + n_rows * (size_t) col_vec[snp];
      unsigned char* dest = reorder ?
        scratch.data() + (size_t) s * n_individuals : col_start + row_offset;
      // with swapped alleles, the REF allele of the file is the alternate
      // allele of the merged locus, so we count it instead (which also gives
      // the right dosage for a genotype with fewer alleles than the ploidy
      // of its individual, e.g. a haploid genotype of a diploid individual)
      int line_status = count_alt_alleles(
        line_start[l], line_len[l], dest, n_individuals, missing,
        swap_vec[snp] ? '0' : '1', chunk_max_alleles);
      if (line_status != GENOTYPES_OK) {
        status = std::max(status, line_status);
        continue;
      }
      if (reorder) {
        for (size_t i = 0; i < n_individuals; i++) {
          col_start[row_vec[i]] = dest[i];
        }
      }
    }
    if (status != GENOTYPES_OK) {
      stop_genotype_status(status, "the VCF " + filename);
    }
    if (chunk_max_alleles > max_ploidy) {
      stop("a genotype of the VCF " + filename + " has more than max_ploidy " +
           "alleles. We estimate max_ploidy from the first variant in the " +
           "vcf file, make sure that variant is representative of ploidy " +
           "(e.g. it is not on a sex chromosome)");
    }
    Rcpp::checkUserInterrupt();
  }
  if (k != n_snps) {
    stop("the VCF " + filename + " has fewer biallelic SNPs than when its " +
         "loci were read (it changed while being merged)");
  }
}
