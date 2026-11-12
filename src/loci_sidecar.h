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
#include <climits>
#include <cctype>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

// A buffered writer of a new sidecar file, which removes the file if it is
// destroyed before finish() (e.g. when an error stops the parsing), as the
// BkWriter of the backing file does. The file is opened in binary mode, so
// that lines end in "\n" on all platforms. With append = true, lines are
// added to the end of an existing file (which is then not removed on error:
// the caller removes it)
class SidecarWriter {
public:
  explicit SidecarWriter(const std::string &filename, bool append = false)
    : filename(filename), finished(false), append(append) {
    fp = std::fopen(filename.c_str(), append ? "ab" : "wb");
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
      if (!append) {
        std::remove(filename.c_str());
      }
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
  bool append;
  std::string buffer;
  void flush() {
    if (!buffer.empty() &&
        (std::fwrite(buffer.data(), 1, buffer.size(), fp) != buffer.size())) {
      Rcpp::stop("error writing the file " + filename);
    }
    buffer.clear();
  }
};

// the length of the prefix "chromosome_", "chr_", "chromosome" or "chr"
// (ignoring case, tried in this order) of a chromosome name, or 0, as
// stripped by cast_chromosome_to_int() in R
inline size_t chromosome_prefix_length(const char* s, size_t len) {
  static const char* prefixes[] = {"chromosome_", "chr_", "chromosome", "chr"};
  for (const char* prefix : prefixes) {
    size_t n = std::strlen(prefix);
    if (len < n) {
      continue;
    }
    size_t i = 0;
    while ((i < n) && (std::tolower((unsigned char) s[i]) == prefix[i])) {
      i++;
    }
    if (i == n) {
      return n;
    }
  }
  return 0;
}

/******************************************************************************
 * The integer codes of chromosomes (chromosomes_as_int)
 *
 * The codes follow cast_chromosome_to_int() in R (used by gt_as_plink(),
 * gt_order_loci() and others), so that the codes are easy to read, with two
 * changes that make them reversible (each name has its own code):
 * - a name that is a positive number, after an optional prefix "chr" or
 *   "chromosome" (see chromosome_prefix_length()), keeps that number as its
 *   code ("7", "chr7", "Chr07"); if several names give the same number
 *   ("1" and "chr1"), only the first to appear keeps it, and the others are
 *   coded as other names; numbers above INT_MAX, and 0, are other names (so
 *   that codes are positive, as the keys of the sidecar);
 * - other names ("chrX", "scaffold_17", "SNP_higher_path_42") are coded after
 *   the largest number, in order of first appearance;
 * - a missing chromosome "." is NA.
 * With millions of chromosomes (e.g. RAD loci assembled de novo), the numbers
 * that are taken are found by sorting pairs of (number, chromosome), rather
 * than with a hash set (about 8 times less memory).
 *
 * @param n the number of chromosomes
 * @param name_of a function name_of(k, start, len) that gives the name of
 *   chromosome k (0-based, in order of first appearance)
 * @return the code of each chromosome (NA_INTEGER for ".")
 */
template <class NameOf>
std::vector<int> chromosome_int_codes(size_t n, NameOf name_of) {
  std::vector<int> codes(n, NA_INTEGER);
  // whether each chromosome is coded by its number (1), missing (2), or an
  // other name (0)
  std::vector<char> kind(n, 0);
  std::vector<std::pair<int, size_t> > numbers;
  for (size_t k = 0; k < n; k++) {
    const char* s;
    size_t len;
    name_of(k, s, len);
    if ((len == 1) && (s[0] == '.')) {
      kind[k] = 2;
      continue;
    }
    size_t i = chromosome_prefix_length(s, len);
    if (i == len) {
      continue;
    }
    int64_t value = 0;
    while ((i < len) && std::isdigit((unsigned char) s[i]) &&
           (value <= INT_MAX)) {
      value = value * 10 + (s[i] - '0');
      i++;
    }
    if ((i == len) && (value >= 1) && (value <= INT_MAX)) {
      numbers.push_back(std::make_pair((int) value, k));
    }
  }
  // the first chromosome (in order of appearance) with each number keeps it
  std::sort(numbers.begin(), numbers.end());
  int64_t max_number = 0;
  for (size_t j = 0; j < numbers.size(); j++) {
    if ((j > 0) && (numbers[j].first == numbers[j - 1].first)) {
      continue;
    }
    codes[numbers[j].second] = numbers[j].first;
    kind[numbers[j].second] = 1;
    max_number = std::max(max_number, (int64_t) numbers[j].first);
  }
  std::vector<std::pair<int, size_t> >().swap(numbers);
  // other names are coded after the largest number
  int64_t next = max_number;
  for (size_t k = 0; k < n; k++) {
    if (kind[k] == 0) {
      next++;
      if (next > INT_MAX) {
        Rcpp::stop("there are too many chromosomes to code them as integers");
      }
      codes[k] = (int) next;
    }
  }
  return codes;
}

// Write the chromosome sidecar (see the top of this file): a line with the
// code and name of each chromosome, in increasing order of code (missing
// chromosomes, coded NA, are left out)
template <class NameOf>
void write_chromosome_sidecar(const std::string &filename,
                              const std::vector<int> &codes,
                              NameOf name_of) {
  std::vector<size_t> order;
  order.reserve(codes.size());
  for (size_t k = 0; k < codes.size(); k++) {
    if (codes[k] != NA_INTEGER) {
      order.push_back(k);
    }
  }
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return codes[a] < codes[b];
  });
  SidecarWriter writer(filename);
  writer.add("int_chromosome\tchromosome");
  writer.end_line();
  for (size_t j = 0; j < order.size(); j++) {
    const char* s;
    size_t len;
    name_of(order[j], s, len);
    writer.add_int(codes[order[j]]);
    writer.tab();
    writer.add(s, len);
    writer.end_line();
  }
  writer.finish();
}

#endif
