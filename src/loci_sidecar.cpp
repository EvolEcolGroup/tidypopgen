// Readers of the sidecar files (see loci_sidecar.h for their format): the
// sparse index of a file, the lines of some keys, all the lines of a file,
// and the copy of the lines of some keys into a new file.
//
// To fetch a few keys from a large file, a sparse index records the key and
// the byte offset of every block-th line (sidecar_index()); a lookup then
// finds the last indexed line at or before the key, seeks there, and reads
// forward. The index is built once per session by R (see R/loci_sidecar.R),
// which keeps it in a cache.

#include "loci_sidecar.h"
#include <vector>
#include <algorithm>
#include <cstring>
#include <climits>

using namespace Rcpp;

// move a file to a byte offset (64-bit, so that files over 2 GB also work on
// Windows, where long has 32 bits)
inline int sidecar_seek(std::FILE* fp, int64_t offset) {
#ifdef _WIN32
  return _fseeki64(fp, (__int64) offset, SEEK_SET);
#else
  return fseeko(fp, (off_t) offset, SEEK_SET);
#endif
}

// A buffered reader of the lines of a sidecar file, which can move to the
// start of any line (given its byte offset)
class SidecarReader {
public:
  explicit SidecarReader(const std::string &filename)
    : filename(filename), buf(1 << 20), start(0), end(0), buf_offset(0),
      at_eof(false) {
    fp = std::fopen(filename.c_str(), "rb");
    if (fp == NULL) {
      stop("could not open the file " + filename);
    }
  }
  ~SidecarReader() {
    if (fp != NULL) {
      std::fclose(fp);
    }
  }
  // move to the line that starts at byte offset
  void seek(int64_t offset) {
    if (sidecar_seek(fp, offset) != 0) {
      stop("error reading the file " + filename);
    }
    start = 0;
    end = 0;
    buf_offset = offset;
    at_eof = false;
  }
  // the next line (without its "\n"; valid until the next call), and the
  // byte offset of its start in the file; false at the end of the file
  bool next_line(const char* &line, size_t &len, int64_t &offset) {
    while (true) {
      const char* nl = static_cast<const char*>(
        std::memchr(buf.data() + start, '\n', end - start));
      if (nl != NULL) {
        line = buf.data() + start;
        len = nl - line;
        offset = buf_offset + start;
        start += len + 1;
        return true;
      }
      if (at_eof) {
        // our files always end with "\n": a last line without it is a sign
        // that the file was cut short (e.g. a full disk, or a copy that was
        // interrupted)
        if (start < end) {
          stop("the file " + filename + " is truncated (its last line does " +
               "not end with a newline)");
        }
        return false;
      }
      // move the incomplete line to the start of the buffer, and read more
      // (growing the buffer if a single line fills it)
      if (start > 0) {
        std::memmove(buf.data(), buf.data() + start, end - start);
        buf_offset += start;
        end -= start;
        start = 0;
      }
      if (end == buf.size()) {
        buf.resize(buf.size() * 2);
      }
      size_t n = std::fread(buf.data() + end, 1, buf.size() - end, fp);
      if (n == 0) {
        if (std::ferror(fp)) {
          stop("error reading the file " + filename);
        }
        at_eof = true;
      }
      end += n;
    }
  }
  // read the first line, and stop if it is not the expected header
  void check_header(const std::string &header) {
    const char* line;
    size_t len;
    int64_t offset;
    if (!next_line(line, len, offset) ||
        (len != header.size()) ||
        (std::memcmp(line, header.data(), len) != 0)) {
      stop("the file " + filename + " is not a sidecar file of tidypopgen " +
           "of the expected type (its first line is not the expected header)");
    }
  }

private:
  std::string filename;
  std::FILE* fp;
  std::vector<char> buf;
  size_t start; // the start of the unread data in buf
  size_t end; // the end of the data in buf
  int64_t buf_offset; // the offset in the file of buf[0]
  bool at_eof;
};

// the integer key of a line (its first field), stopping if it is not a
// positive integer
inline int parse_key(const char* line, size_t len, const std::string &path) {
  size_t i = 0;
  int64_t key = 0;
  while ((i < len) && (line[i] >= '0') && (line[i] <= '9') &&
         (key <= INT_MAX)) {
    key = key * 10 + (line[i] - '0');
    i++;
  }
  if ((i == 0) || (key < 1) || (key > INT_MAX) ||
      ((i < len) && (line[i] != '\t'))) {
    stop("the file " + path + " has a line whose first field is not a " +
         "positive integer key");
  }
  return (int) key;
}

// an integer field (of length len), or NA for a missing value "."; it stops
// if the field is not an integer
inline int parse_int_field(const char* p, size_t len, const std::string &path) {
  if ((len == 0) || ((len == 1) && (p[0] == '.'))) {
    return NA_INTEGER;
  }
  size_t i = 0;
  bool negative = (p[0] == '-');
  if (negative) {
    i++;
  }
  int64_t x = 0;
  size_t first_digit = i;
  while ((i < len) && (p[i] >= '0') && (p[i] <= '9') && (x <= INT_MAX)) {
    x = x * 10 + (p[i] - '0');
    i++;
  }
  if ((i == first_digit) || (i != len) || (x > INT_MAX)) {
    stop("the file " + path + " has a value that is not an integer");
  }
  return (int) (negative ? -x : x);
}

// split a line into its first n_fields fields (1-based: field 1 is the key);
// missing fields (a line with fewer fields) get a NULL start
inline void split_fields(const char* line, size_t len, size_t n_fields,
                         std::vector<const char*> &starts,
                         std::vector<size_t> &lens) {
  const char* p = line;
  const char* line_end = line + len;
  for (size_t f = 1; f <= n_fields; f++) {
    if (p == NULL) {
      starts[f] = NULL;
      lens[f] = 0;
      continue;
    }
    const char* tab = static_cast<const char*>(
      std::memchr(p, '\t', line_end - p));
    starts[f] = p;
    lens[f] = (tab == NULL) ? (size_t) (line_end - p) : (size_t) (tab - p);
    p = (tab == NULL) ? NULL : tab + 1;
  }
}

// the columns requested from a sidecar file (fields, 1-based, with their
// types), to be filled for n items
struct SidecarColumns {
  std::vector<int> fields;
  std::vector<bool> is_int;
  std::vector<SEXP> cols;
  List out;
  size_t max_field;
  SidecarColumns(IntegerVector fields_r, CharacterVector types, size_t n)
    : out(fields_r.size()), max_field(1) {
    if (fields_r.size() != types.size()) {
      stop("fields and types must have the same length");
    }
    for (R_xlen_t j = 0; j < fields_r.size(); j++) {
      if ((fields_r[j] == NA_INTEGER) || (fields_r[j] < 1)) {
        stop("fields must be positive integers");
      }
      fields.push_back(fields_r[j]);
      max_field = std::max(max_field, (size_t) fields_r[j]);
      std::string type = as<std::string>(types[j]);
      if ((type != "integer") && (type != "character")) {
        stop("types must be 'integer' or 'character'");
      }
      is_int.push_back(type == "integer");
      if (is_int[j]) {
        out[j] = IntegerVector(n, NA_INTEGER);
      } else {
        CharacterVector col(n);
        for (size_t i = 0; i < n; i++) {
          SET_STRING_ELT(col, i, NA_STRING);
        }
        out[j] = col;
      }
      // (the columns are protected by out)
      cols.push_back(out[j]);
    }
  }
  // store the fields of a line as item i
  void set(size_t i, const std::vector<const char*> &starts,
           const std::vector<size_t> &lens, const std::string &path) {
    for (size_t j = 0; j < fields.size(); j++) {
      const char* p = starts[fields[j]];
      size_t len = lens[fields[j]];
      if (is_int[j]) {
        INTEGER(cols[j])[i] = (p == NULL) ? NA_INTEGER :
          parse_int_field(p, len, path);
      } else if ((p != NULL) && !((len == 1) && (p[0] == '.'))) {
        SET_STRING_ELT(cols[j], i, Rf_mkCharLen(p, len));
      }
    }
  }
  // copy item i_from to item i
  void copy(size_t i, size_t i_from) {
    for (size_t j = 0; j < fields.size(); j++) {
      if (is_int[j]) {
        INTEGER(cols[j])[i] = INTEGER(cols[j])[i_from];
      } else {
        SET_STRING_ELT(cols[j], i, STRING_ELT(cols[j], i_from));
      }
    }
  }
};

/******************************************************************************
 * Build the sparse index of a sidecar file
 *
 * It reads the file once, checking its header and that its keys are positive
 * integers that strictly increase, and records the key and the byte offset
 * of every block-th line (the first line of each block of lines).
 *
 * @param path the file
 * @param header its expected first line
 * @param block the number of lines per block
 * @return a list with the `keys` and `offsets` (as doubles, which hold any
 *   offset of a file below 8 PB exactly) of the first line of each block,
 *   and `n`, the number of lines (without the header)
 */
// [[Rcpp::export]]
List sidecar_index(std::string path, std::string header, int block) {
  if (block < 1) {
    stop("block must be a positive integer");
  }
  SidecarReader reader(path);
  reader.check_header(header);
  std::vector<int> keys;
  std::vector<double> offsets;
  const char* line;
  size_t len;
  int64_t offset;
  int64_t n = 0;
  int previous = 0;
  while (reader.next_line(line, len, offset)) {
    int key = parse_key(line, len, path);
    if (key <= previous) {
      stop("the keys of the file " + path + " are not strictly increasing " +
           "(line " + std::to_string(n + 2) + ")");
    }
    if (n % block == 0) {
      keys.push_back(key);
      offsets.push_back((double) offset);
    }
    previous = key;
    n++;
    if ((n & 0xFFFFF) == 0) {
      checkUserInterrupt();
    }
  }
  return List::create(_["keys"] = IntegerVector(keys.begin(), keys.end()),
                      _["offsets"] = NumericVector(offsets.begin(),
                                                   offsets.end()),
                      _["n"] = (double) n);
}

/******************************************************************************
 * Fetch the lines of some keys from a sidecar file
 *
 * The keys can be in any order, with repeats and NAs. They are looked up in
 * increasing order, so that each block of the file is read at most once: for
 * each key, if the reader is before the block of the key (the last indexed
 * line at or before it), it seeks to the block, then reads forward to the
 * key. Keys in increasing order (e.g. all the loci of a gen_tibble read from
 * a VCF) are read in a single pass through the file.
 *
 * Optionally, an integer field is compared with expected values (the
 * positions of the loci in the gen_tibble, to check that the file belongs to
 * it): NAs are equal to each other.
 *
 * @param path the file
 * @param keys the keys to fetch
 * @param index_keys,index_offsets the sparse index (from sidecar_index())
 * @param fields the fields to return (1-based; field 1 is the key)
 * @param types the type of each field ("character" or "integer")
 * @param check_field the field to compare with check_values (0 for none)
 * @param check_values the expected value of check_field for each key
 * @return a list with a vector for each field (NA for keys that are not
 *   found, and for missing values "."), `found` (whether each key was
 *   found), `n_mismatch` (the number of keys whose check_field differs from
 *   check_values), and `first_mismatch` (the index, 1-based, of such a key,
 *   or NA) with the value of its check_field in the file
 *   (`first_mismatch_value`)
 */
// [[Rcpp::export]]
List sidecar_fetch(std::string path, IntegerVector keys,
                   IntegerVector index_keys, NumericVector index_offsets,
                   IntegerVector fields, CharacterVector types,
                   int check_field, IntegerVector check_values) {
  size_t n = keys.size();
  SidecarColumns columns(fields, types, n);
  LogicalVector found(n, false);
  bool check = (check_field > 0);
  if (check && ((size_t) check_values.size() != n)) {
    stop("check_values must have one value per key");
  }
  size_t n_split = std::max(columns.max_field, (size_t) check_field);
  std::vector<const char*> starts(n_split + 1);
  std::vector<size_t> lens(n_split + 1);

  // the keys to look up (not NA), in increasing order (sorting is skipped if
  // they already are)
  std::vector<int> order;
  order.reserve(n);
  bool sorted = true;
  for (size_t i = 0; i < n; i++) {
    if (keys[i] != NA_INTEGER) {
      if (!order.empty() && (keys[i] < keys[order.back()])) {
        sorted = false;
      }
      order.push_back(i);
    }
  }
  if (!sorted) {
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
      return keys[a] < keys[b];
    });
  }

  SidecarReader reader(path);
  // the current line of the reader (the first line with a key at or after
  // the previous key looked up), and whether there is one (false before the
  // first read, and at the end of the file)
  const char* line = NULL;
  size_t len = 0;
  int64_t offset;
  int line_key = 0;
  bool started = false;
  bool at_end = false;
  auto read_line = [&]() {
    if (reader.next_line(line, len, offset)) {
      line_key = parse_key(line, len, path);
    } else {
      at_end = true;
    }
  };

  double n_mismatch = 0;
  int first_mismatch = NA_INTEGER;
  int first_mismatch_value = NA_INTEGER;
  int previous_key = NA_INTEGER; // the previous key looked up
  int previous_i = -1; // its index
  int previous_value = NA_INTEGER; // the value of its check_field
  for (size_t k = 0; k < order.size(); k++) {
    size_t i = order[k];
    int key = keys[i];
    int value = NA_INTEGER; // the value of check_field in the file
    if ((previous_i >= 0) && (key == previous_key)) {
      // a repeated key: copy what we found for it
      if (found[previous_i]) {
        columns.copy(i, previous_i);
        found[i] = true;
        value = previous_value;
      }
    } else {
      // the block of the key: the last indexed line at or before it
      IntegerVector::iterator block = std::upper_bound(
        index_keys.begin(), index_keys.end(), key);
      if (block != index_keys.begin()) {
        size_t b = (block - index_keys.begin()) - 1;
        // seek to the block if the reader is before it
        if (!started || (!at_end && (line_key < index_keys[b]))) {
          reader.seek((int64_t) index_offsets[b]);
          at_end = false;
          started = true;
          read_line();
        }
        while (!at_end && (line_key < key)) {
          read_line();
        }
        if (!at_end && (line_key == key)) {
          split_fields(line, len, n_split, starts, lens);
          columns.set(i, starts, lens, path);
          found[i] = true;
          if (check) {
            value = (starts[check_field] == NULL) ? NA_INTEGER :
              parse_int_field(starts[check_field], lens[check_field], path);
          }
        }
      }
    }
    if (check && found[i] && (value != check_values[i])) {
      n_mismatch++;
      if ((first_mismatch == NA_INTEGER) || ((int) i + 1 < first_mismatch)) {
        first_mismatch = i + 1;
        first_mismatch_value = value;
      }
    }
    previous_key = key;
    previous_i = i;
    previous_value = value;
    if ((k & 0xFFFFF) == 0) {
      checkUserInterrupt();
    }
  }
  List res = columns.out;
  res.push_back(found, "found");
  res.push_back(n_mismatch, "n_mismatch");
  res.push_back(first_mismatch, "first_mismatch");
  res.push_back(first_mismatch_value, "first_mismatch_value");
  return res;
}

/******************************************************************************
 * Read all the lines of a sidecar file
 *
 * @param path the file
 * @param header its expected first line
 * @param fields the fields to return (1-based; field 1 is the key)
 * @param types the type of each field ("character" or "integer")
 * @param n the number of lines (without the header), e.g. from
 *   sidecar_index(), or -1 to count them first (with an extra pass through
 *   the file)
 * @return a list with a vector for each field
 */
// [[Rcpp::export]]
List sidecar_read(std::string path, std::string header, IntegerVector fields,
                  CharacterVector types, double n = -1) {
  const char* line;
  size_t len;
  int64_t offset;
  if (n < 0) {
    SidecarReader counter(path);
    counter.check_header(header);
    n = 0;
    while (counter.next_line(line, len, offset)) {
      n++;
    }
  }
  SidecarColumns columns(fields, types, (size_t) n);
  std::vector<const char*> starts(columns.max_field + 1);
  std::vector<size_t> lens(columns.max_field + 1);
  SidecarReader reader(path);
  reader.check_header(header);
  size_t i = 0;
  while (reader.next_line(line, len, offset)) {
    if (i >= (size_t) n) {
      stop("the file " + path + " has more lines than expected");
    }
    split_fields(line, len, columns.max_field, starts, lens);
    columns.set(i, starts, lens, path);
    i++;
    if ((i & 0xFFFFF) == 0) {
      checkUserInterrupt();
    }
  }
  if (i != (size_t) n) {
    stop("the file " + path + " has fewer lines than expected");
  }
  return columns.out;
}

/******************************************************************************
 * Copy the lines of some keys of a sidecar file into a new file
 *
 * This writes the sidecar of a new backing file (e.g. from
 * gt_update_backingfile() or rbind()), with the lines of the loci that it
 * keeps. The keys to keep are held in a bitset (one bit per possible key, so
 * N/8 bytes), and the file is streamed once, in the order of its keys
 * (whatever the order of the loci in the gen_tibble). The new file is removed
 * if an error stops the copy.
 *
 * @param path_in the sidecar to copy from
 * @param path_out the new sidecar
 * @param header the expected first line of path_in (copied to path_out)
 * @param keys the keys to keep (positive integers, in any order)
 * @return a list with `n_keys`, the number of distinct keys, and
 *   `n_written`, the number of lines written (fewer than n_keys if some keys
 *   are not in path_in)
 */
// [[Rcpp::export]]
List sidecar_subset(std::string path_in, std::string path_out,
                    std::string header, IntegerVector keys) {
  int max_key = 0;
  for (R_xlen_t i = 0; i < keys.size(); i++) {
    if ((keys[i] == NA_INTEGER) || (keys[i] < 1)) {
      stop("keys must be positive integers");
    }
    max_key = std::max(max_key, keys[i]);
  }
  std::vector<uint64_t> keep(((size_t) max_key >> 6) + 1, 0);
  double n_keys = 0;
  for (R_xlen_t i = 0; i < keys.size(); i++) {
    uint64_t bit = (uint64_t) 1 << (keys[i] & 63);
    uint64_t &word = keep[keys[i] >> 6];
    if (!(word & bit)) {
      word |= bit;
      n_keys++;
    }
  }
  SidecarReader reader(path_in);
  reader.check_header(header);
  SidecarWriter writer(path_out);
  writer.add(header);
  writer.end_line();
  const char* line;
  size_t len;
  int64_t offset;
  double n_written = 0;
  int previous = 0;
  while (reader.next_line(line, len, offset)) {
    int key = parse_key(line, len, path_in);
    if (key <= previous) {
      stop("the keys of the file " + path_in + " are not strictly " +
           "increasing");
    }
    previous = key;
    // the keys increase, so there is nothing more to copy
    if (key > max_key) {
      break;
    }
    if ((keep[key >> 6] >> (key & 63)) & 1) {
      writer.add(line, len);
      writer.end_line();
      n_written++;
    }
  }
  writer.finish();
  return List::create(_["n_keys"] = n_keys, _["n_written"] = n_written);
}
