#include <Rcpp.h>
#include <cstring>
#include <cstdint>
#include <vector>
using namespace Rcpp;

/******************************************************************************/
// Membership of the elements of a character vector in a small table of
// strings, as with %in%, used to check and harmonise the alleles of the loci
// table (see check_allele_alphabet() and harmonise_missing_values()).
//
// With tens of millions of loci, %in% (or vctrs::vec_in()) allocates several
// vectors as long as the loci table, which add up to more than a GB of
// transient memory. Here each element is checked in a single pass, without
// allocating anything per element.
//
// R keeps a single copy (CHARSXP) of each string with a given encoding, and
// a single copy of each ASCII string (whatever its declared encoding), so an
// element is in the table if it is the same CHARSXP as one of the strings in
// the table. Otherwise (e.g. for the same non-ASCII string in two encodings)
// we compare the characters, translated to UTF-8 (or as bytes, if either
// string is marked as "bytes").
//
// As the alleles of a loci table take few distinct values, the result for
// each distinct CHARSXP is computed once and kept in a small hash table keyed
// by its address, so that each element then costs a single lookup (scanning
// the strings of the table instead is several times slower, as the branches
// depend on the alleles, which are in random order).

class StringTable {
public:
  StringTable(const CharacterVector& table) : n_cached(0) {
    for (R_xlen_t j = 0; j < table.size(); j++) {
      strings.push_back(table[j]);
    }
    for (size_t h = 0; h < cache_size; h++) {
      keys[h] = NULL;
      in_table[h] = 0;
    }
  }

  // is s (a CHARSXP) in the table?
  inline bool contains(SEXP s) {
    // find s in the hash table (open addressing with linear probing)
    size_t h = hash(s);
    while (keys[h] != NULL) {
      if (keys[h] == s) return in_table[h];
      h = (h + 1) & (cache_size - 1);
    }
    bool res = compare(s);
    // keep the result, leaving at least half of the hash table empty (if
    // there are more distinct strings, the others are compared every time)
    if (n_cached < cache_size / 2) {
      keys[h] = s;
      in_table[h] = res;
      n_cached++;
    }
    return res;
  }

private:
  static const size_t cache_size = 256; // a power of 2
  std::vector<SEXP> strings;
  SEXP keys[cache_size];
  char in_table[cache_size];
  size_t n_cached;

  // multiplicative (Fibonacci) hash of the address of a CHARSXP
  inline size_t hash(SEXP s) {
    uint64_t a = (uint64_t)(uintptr_t) s;
    return (size_t)((a * 0x9E3779B97F4A7C15ULL) >> 56) & (cache_size - 1);
  }

  // is s the same string as one in the table?
  bool compare(SEXP s) {
    for (size_t j = 0; j < strings.size(); j++) {
      if (s == strings[j]) return true;
    }
    // NA is a single CHARSXP, so it only matches itself
    if (s == NA_STRING) return false;
    for (size_t j = 0; j < strings.size(); j++) {
      SEXP t = strings[j];
      if (t == NA_STRING) continue;
      if (Rf_getCharCE(s) == CE_BYTES || Rf_getCharCE(t) == CE_BYTES) {
        if (std::strcmp(CHAR(s), CHAR(t)) == 0) return true;
      } else if (std::strcmp(Rf_translateCharUTF8(s),
                             Rf_translateCharUTF8(t)) == 0) {
        return true;
      }
    }
    return false;
  }
};

/******************************************************************************/
// all(x %in% table), stopping at the first element of x not in table

// [[Rcpp::export]]
bool all_in_cpp(const CharacterVector& x, const CharacterVector& table) {
  StringTable lookup(table);
  // read the elements directly (STRING_ELT() is a function call per element)
  const SEXP* p_x = STRING_PTR_RO(x);
  R_xlen_t n = x.size();
  for (R_xlen_t i = 0; i < n; i++) {
    if (!lookup.contains(p_x[i])) return false;
  }
  return true;
}

/******************************************************************************/
// which(x %in% table), as 1-based indices (doubles if x is a long vector)

// [[Rcpp::export]]
SEXP which_in_cpp(const CharacterVector& x, const CharacterVector& table) {
  StringTable lookup(table);
  std::vector<R_xlen_t> found;
  const SEXP* p_x = STRING_PTR_RO(x);
  R_xlen_t n = x.size();
  for (R_xlen_t i = 0; i < n; i++) {
    if (lookup.contains(p_x[i])) found.push_back(i + 1);
  }
  if (x.size() <= INT_MAX) {
    return IntegerVector(found.begin(), found.end());
  } else {
    return NumericVector(found.begin(), found.end());
  }
}

/******************************************************************************/
// The index (1-based) of the first element of a character vector that is a
// duplicate of an earlier one, or 0 if there is none, as anyDuplicated(), or
// -1 if this can not be decided here (and anyDuplicated() should be used),
// to check that the names of the loci are unique (see any_duplicated_names()).
//
// anyDuplicated() hashes the addresses of the CHARSXPs, with a table of 2^k
// ints for at least 2 * n elements (512 MB for 40M names). Here the table is
// up to 70% full (256 MB for 40M names), and the slot of each element is
// prefetched a few elements ahead, which makes it about 3 times faster. The
// addresses identify the strings if all the elements are in the native
// encoding (as ASCII strings always are, as they are never marked): R keeps
// a single CHARSXP for each string in each encoding, but the same text can
// be in two encodings (e.g. UTF-8 and latin1), so other vectors are left to
// anyDuplicated() (which translates them).

namespace {
inline void prefetch_slot(const void* a) {
#if defined(__GNUC__) || defined(__clang__)
  __builtin_prefetch(a);
#else
  (void) a;
#endif
}

inline uint64_t string_address_hash(SEXP s) {
  uint64_t a = (uint64_t)(uintptr_t) s;
  a ^= a >> 33;
  a *= 0xff51afd7ed558ccdULL;
  a ^= a >> 33;
  return a;
}
} // namespace

// [[Rcpp::export]]
double any_duplicated_strings_cpp(const CharacterVector& x) {
  const SEXP* p = STRING_PTR_RO(x);
  R_xlen_t n = x.size();
  if ((uint64_t) n >= (uint64_t) UINT32_MAX) return -1;
  for (R_xlen_t i = 0; i < n; i++) {
    if (p[i] != NA_STRING && Rf_getCharCE(p[i]) != CE_NATIVE) return -1;
  }
  // a table of 2^k slots, at most 70% full, holding the index + 1 of each
  // element (0 for an empty slot), with linear probing; the slot is given by
  // the top k bits of the hash
  int k = 4;
  while ((double)((uint64_t) 1 << k) * 0.7 < (double) n) k++;
  const uint64_t mask = ((uint64_t) 1 << k) - 1;
  const int shift = 64 - k;
  std::vector<uint32_t> slots((size_t) 1 << k, 0);
  const R_xlen_t ahead = 16;
  for (R_xlen_t i = 0; i < n; i++) {
    if (i + ahead < n) {
      prefetch_slot(&slots[string_address_hash(p[i + ahead]) >> shift]);
    }
    SEXP s = p[i];
    uint64_t h = string_address_hash(s) >> shift;
    while (slots[h] != 0) {
      if (p[slots[h] - 1] == s) return (double)(i + 1);
      h = (h + 1) & mask;
    }
    slots[h] = (uint32_t)(i + 1);
  }
  return 0;
}

/******************************************************************************/
