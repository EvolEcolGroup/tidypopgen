#ifndef TIDYPOPGEN_LOCI_SIDECAR_H
#define TIDYPOPGEN_LOCI_SIDECAR_H

// Sidecar files hold the original names of the loci (and chromosomes) that a
// gen_tibble stores as integers (names_as_int, chromosomes_as_int). They are
// plain tab separated text: a header line, then one line per locus (or
// chromosome) whose first field is its integer key, with the keys strictly
// increasing; there is no quoting, "." is a missing value, and lines end in
// "\n" on all platforms. They sit next to the backing file (foo.bk ->
// foo.loci.tsv, foo.chromosomes.tsv), and are never edited: functions that
// write a new backing file write new sidecars for it (with sidecar_subset()).
// The readers are in loci_sidecar.cpp; this header has the writer, which the
// VCF parser uses to stream the locus sidecar while it reads the file.

#include <Rcpp.h>
#include <cstdio>
#include <cstdint>
#include <string>

// A buffered writer of a new sidecar file, which removes the file if it is
// destroyed before finish() (e.g. when an error stops the parsing), as the
// BkWriter of the backing file does. The file is opened in binary mode, so
// that lines end in "\n" on all platforms.
class SidecarWriter {
public:
  explicit SidecarWriter(const std::string &filename)
    : filename(filename), finished(false) {
    fp = std::fopen(filename.c_str(), "wb");
    if (fp == NULL) {
      Rcpp::stop("could not create the file " + filename);
    }
    // (a little more than the flush size, so that the last line of a buffer
    // rarely needs a reallocation)
    buffer.reserve(flush_size + 4096);
  }
  ~SidecarWriter() {
    if (!finished) {
      if (fp != NULL) {
        std::fclose(fp);
      }
      std::remove(filename.c_str());
    }
  }
  void add(const char* x, size_t len) { buffer.append(x, len); }
  void add(const std::string &x) { buffer.append(x); }
  // append an integer in decimal (as as.character() writes an integer)
  void add_int(int64_t x) {
    char digits[24];
    int n = 0;
    uint64_t u = (x < 0) ? (uint64_t) (-(x + 1)) + 1 : (uint64_t) x;
    do {
      digits[n++] = '0' + (char) (u % 10);
      u /= 10;
    } while (u > 0);
    if (x < 0) {
      buffer.push_back('-');
    }
    while (n > 0) {
      buffer.push_back(digits[--n]);
    }
  }
  void tab() { buffer.push_back('\t'); }
  // end a line, writing the buffer to the file once it is large enough
  void end_line() {
    buffer.push_back('\n');
    if (buffer.size() >= flush_size) {
      flush();
    }
  }
  void finish() {
    flush();
    if (std::fclose(fp) != 0) {
      fp = NULL;
      Rcpp::stop("error writing the file " + filename);
    }
    fp = NULL;
    finished = true;
  }

private:
  static const size_t flush_size = 1 << 22;
  std::string filename;
  std::FILE* fp;
  bool finished;
  std::string buffer;
  void flush() {
    if (!buffer.empty() &&
        (std::fwrite(buffer.data(), 1, buffer.size(), fp) != buffer.size())) {
      Rcpp::stop("error writing the file " + filename);
    }
    buffer.clear();
  }
};

#endif
