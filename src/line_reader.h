#ifndef TIDYPOPGEN_LINE_READER_H
#define TIDYPOPGEN_LINE_READER_H

// The reader of the lines of a text file (plain, gzipped, or bgzipped and
// decompressed in parallel) used by the VCF parser (vcf_parser_fast.cpp) and
// by the PLINK .ped reader (ped_to_fbm.cpp), a chunk of whole lines at a time,
// and a buffer of strings (StringBuffer) used by both to collect names

#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdio>
#include <zlib.h>

using Rcpp::stop;

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
  Rcpp::CharacterVector to_r() {
    Rcpp::CharacterVector res(size());
    for (size_t i = 0; i < size(); i++) {
      SET_STRING_ELT(res, i, Rf_mkCharLen(chars.data() + starts[i],
                                          starts[i + 1] - starts[i]));
    }
    std::string().swap(chars);
    std::vector<size_t>().swap(starts);
    return res;
  }
};

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

// Function to check whether a file starts with the two magic bytes of gzip
// (31, 139), as R's gzfile() does to decide whether to decompress a file
inline bool has_gz_magic(const std::string &filename) {
  std::FILE* f = std::fopen(filename.c_str(), "rb");
  if (f == NULL) {
    return false;
  }
  unsigned char magic[2];
  bool is_gz = (std::fread(magic, 1, 2, f) == 2) && (magic[0] == 31) &&
    (magic[1] == 139);
  std::fclose(f);
  return is_gz;
}

// A reader of the lines of a VCF (or of another text file, such as a PLINK
// .ped file; plain or gzipped), a chunk at a time: each
// chunk holds whole lines, each ended by '\0' (written in place of the
// newline, and of a carriage return before it), so that they can be parsed
// as C strings. The lines of a chunk are valid until the next chunk is read.
// Files compressed with bgzip (BGZF, made of independent blocks of at most
// 64 KB) are read a batch of blocks at a time, and the blocks of a batch are
// decompressed in parallel (with n_threads threads); other gzipped files are
// decompressed sequentially with zlib.
class VcfLineReader {
public:
  // chunk_bytes is the size of a chunk (larger chunks hold more lines to
  // parse in parallel); with sniff_gz = true, a file is read as gzipped if
  // it starts with the gzip magic bytes, whatever its extension (as R's
  // gzfile() does), otherwise only if its name ends in .gz
  VcfLineReader(const std::string &filename, int n_threads,
                size_t chunk_bytes = 1 << 23, bool sniff_gz = false)
    : fp(NULL), gz(NULL), is_bgzf(false), n_threads(n_threads),
      chunk_bytes(chunk_bytes), buffer(chunk_bytes), carry_start(0),
      n_carry(0), at_eof(false) {
    // (a thread count that is not positive, e.g. NA, would make OpenMP
    // fail, or try to start a huge number of threads)
    if (n_threads < 1) {
      stop("'n_threads' should be a positive integer");
    }
    if (is_gz_file(filename) || (sniff_gz && has_gz_magic(filename))) {
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
      // grow the buffer once it is full of a line longer than it (a read can
      // fill only part of the buffer, see read(); a batch of BGZF blocks
      // grows the buffer as needed)
      if (!is_bgzf && (n_data == buffer.size())) {
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
  std::FILE* fp;
  gzFile gz;
  bool is_bgzf;
  int n_threads;
  size_t chunk_bytes;
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
      // gzread() reads at most INT_MAX bytes at a time (larger requests are
      // an error, and its length is an unsigned int), so a large buffer
      // (e.g. the chunks of a .ped file with millions of loci, read on many
      // threads) is filled 1 GB at a time; next_chunk() reads again until
      // it finds a newline
      const size_t max_gz_read = (size_t) 1 << 30;
      unsigned int n_request = (unsigned int) std::min(n, max_gz_read);
      int n_read = gzread(gz, dest, n_request);
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

#endif // TIDYPOPGEN_LINE_READER_H
