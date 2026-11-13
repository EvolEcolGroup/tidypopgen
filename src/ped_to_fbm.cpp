#include <Rcpp.h>
#include <vector>
#include <string>
#include <cstring>
#include <cstdint>
#include <unordered_map>
#include <bigstatsr/BMCodeAcc.h>
#include "line_reader.h"

using namespace Rcpp;

/******************************************************************************
 * Read a PLINK .ped file straight into the backing file of an FBM
 *
 * A .ped file has a line per individual: six fields as in a .fam file
 * (family, individual, paternal and maternal IDs, sex and phenotype), then
 * the two alleles of each locus (in the order of the .map file). The FBM
 * stores, for each individual and locus, the count of the second allele of
 * the locus (0, 1 or 2; 3 for a missing genotype), the genotypes of a locus
 * one after the other (column major).
 *
 * The coding of each locus follows snpStats::read.pedfile(), which the
 * previous R reader (read_pedfile()) was modified from: the first allele
 * (allele_ref) of a locus is the first non-missing allele found in the file
 * (reading the individuals in order, and the two alleles of each in order),
 * and the second allele (allele_alt) the first one that differs from it. The
 * coding of a genotype is therefore final as soon as it is read: it only
 * depends on the genotypes before it. The only exception are multi-allelic
 * loci, whose genotypes are all set to missing at the end: their columns are
 * rewritten once the whole file has been read. So the genotypes are written
 * straight into the FBM in a single pass over the file (after a quick pass to
 * count the individuals, which the FBM needs to be created), and the memory
 * used does not grow with the number of individuals times loci.
 *
 * The genotype of an individual is coded from flags `one` (an allele equals
 * the first allele of the locus) and `two` (an allele equals the second):
 * one only gives 0, both 1, two only 2, and neither (both alleles missing)
 * is missing. As in snpStats, a genotype with a single missing allele counts
 * its other allele (e.g. "A 0" is coded as "A A"; we keep this quirk so
 * that the results do not change). Any third allele, in either position of a
 * genotype (e.g. "C A" or "A C", with alleles A and G), makes the locus
 * multi-allelic (snpStats, and the previous reader, missed it in the second
 * position if the first allele matched, and coded "A C" as "A A").
 *
 * The lines are read a chunk at a time (see VcfLineReader in line_reader.h).
 * The lines of a chunk are split in parallel (with n_threads threads), each
 * allele being turned into a small integer id; the genotypes of the chunk
 * are then coded in parallel over loci (the coding of a locus depends on the
 * individuals before, but not on other loci), each thread writing the
 * genotypes of its loci straight into the (memory mapped) FBM.
 *****************************************************************************/

// whether a character is white space for the purpose of finding blank lines
// (as "^\\s*$" in R)
inline bool is_space_char(char c) {
  return (c == ' ') || (c == '\t') || (c == '\r') || (c == '\n') ||
    (c == '\f') || (c == '\v');
}

// whether a line only contains white space (blank lines are skipped)
inline bool is_blank_line(const char* line, size_t len) {
  for (size_t i = 0; i < len; i++) {
    if (!is_space_char(line[i])) {
      return false;
    }
  }
  return true;
}

// whether a character separates the fields of a line (runs of spaces and
// tabs, as PLINK reads .ped files)
inline bool is_field_sep(char c) {
  return (c == ' ') || (c == '\t');
}

/******************************************************************************
 * Count the individuals (the lines that are not blank) of a .ped file
 *
 * @param filename the .ped file (plain or gzipped)
 * @param n_threads the number of threads (to decompress a bgzipped file)
 * @return the number of lines that are not blank
 *****************************************************************************/
// [[Rcpp::export]]
double count_ped_individuals(std::string filename, int n_threads = 1) {
  // (gzipped files are recognised from their first bytes, as by gzfile())
  VcfLineReader reader(filename, n_threads, 1 << 23, true);
  std::vector<char*> line_start;
  std::vector<size_t> line_len;
  double n = 0;
  while (reader.next_chunk(line_start, line_len)) {
    // (reading a large file takes a while: allow the user to interrupt it,
    // as in ped_to_fbm(); the reader closes the file as the exception
    // unwinds the stack)
    checkUserInterrupt();
    for (size_t l = 0; l < line_start.size(); l++) {
      if (!is_blank_line(line_start[l], line_len[l])) {
        n++;
      }
    }
  }
  return n;
}

// The table of the ids of the alleles. An allele of a single character c is
// its own id (1 to 255, so that nothing has to be looked up), and longer
// alleles are numbered from 256 (looked up in a hash table; they are rare,
// as only alleles listed in valid_alleles are accepted later). Missing
// alleles (na_strings) have id 0.
class AlleleIds {
public:
  explicit AlleleIds(const CharacterVector &na_strings) {
    for (int c = 0; c < 256; c++) {
      single[c] = (uint16_t) c;
    }
    // (a NUL character can not be an allele: it is missing)
    single[0] = 0;
    for (R_xlen_t k = 0; k < na_strings.size(); k++) {
      if (STRING_ELT(na_strings, k) == NA_STRING) {
        continue;
      }
      std::string na(na_strings[k]);
      if (na.size() == 1) {
        single[(unsigned char) na[0]] = 0;
      } else if (na.size() > 1) {
        multi[na] = 0;
      }
    }
  }

  // the id of an allele (of len characters starting at p); it returns
  // too_many if there are more distinct long alleles than ids. Long alleles
  // are added to the table under an OpenMP critical section, so that this
  // can be called by several threads
  uint16_t id(const char* p, size_t len) {
    if (len == 1) {
      return single[(unsigned char) p[0]];
    }
    uint16_t res;
#pragma omp critical(ped_allele_ids)
{
    std::string allele(p, len);
    std::unordered_map<std::string, uint16_t>::iterator it =
      multi.find(allele);
    if (it != multi.end()) {
      res = it->second;
    } else if (names.size() >= max_multi) {
      res = too_many;
    } else {
      res = (uint16_t) (256 + names.size());
      multi[allele] = res;
      names.push_back(allele);
    }
}
    return res;
  }

  // the allele with a given id (NA for 0, i.e. no allele)
  SEXP name(uint16_t id) const {
    if (id == 0) {
      return NA_STRING;
    }
    if (id < 256) {
      char c = (char) id;
      return Rf_mkCharLen(&c, 1);
    }
    const std::string &allele = names[id - 256];
    return Rf_mkCharLen(allele.data(), allele.size());
  }

  // the id returned when there are too many distinct long alleles
  static const uint16_t too_many = 65535;

private:
  static const size_t max_multi = 65535 - 256;
  uint16_t single[256];
  std::unordered_map<std::string, uint16_t> multi;
  std::vector<std::string> names;
};

// What we find on a line of a .ped file: the start and length of its first
// six fields, its number of fields, and the first locus whose allele could
// not be given an id (-1 if none)
struct PedLine {
  const char* fam_start[6];
  size_t fam_len[6];
  size_t n_fields;
  long bad_locus;
};

// Split a line of a .ped file: store its first six fields in res, and the ids
// of the 2 * n_loci alleles in ids (two per locus, in order); fields beyond
// these are only counted (a line with too few or too many fields is an error,
// reported by the caller)
inline void split_ped_line(const char* p, size_t len, size_t n_loci,
                           AlleleIds &allele_ids, uint16_t* ids,
                           PedLine &res) {
  const char* end = p + len;
  const size_t n_wanted = 6 + 2 * n_loci;
  size_t field = 0;
  res.bad_locus = -1;
  while (true) {
    // skip the separators before the field
    while ((p < end) && is_field_sep(*p)) {
      p++;
    }
    if (p == end) {
      break;
    }
    const char* start = p;
    while ((p < end) && !is_field_sep(*p)) {
      p++;
    }
    size_t field_len = p - start;
    if (field < 6) {
      res.fam_start[field] = start;
      res.fam_len[field] = field_len;
    } else if (field < n_wanted) {
      uint16_t id = allele_ids.id(start, field_len);
      if ((id == AlleleIds::too_many) && (res.bad_locus < 0)) {
        res.bad_locus = (field - 6) / 2;
      }
      ids[field - 6] = id;
    }
    field++;
  }
  res.n_fields = field;
}

// the label of an individual in error messages: its ID (the second field of
// its line, if any) and the line of the file it is on
inline std::string individual_label(const PedLine &line, double line_number) {
  std::string label = "the individual ";
  if (line.n_fields >= 2) {
    label += std::string(line.fam_start[1], line.fam_len[1]) + " ";
  }
  char buf[64];
  std::snprintf(buf, sizeof(buf), "(line %.0f of the ped file)", line_number);
  return label + buf;
}

/******************************************************************************
 * Write the genotypes of a .ped file into an FBM
 *
 * @param filename the .ped file (plain or gzipped)
 * @param BM the FBM.code256 to fill, with a row per individual (the lines
 *   of the file that are not blank, see count_ped_individuals()) and a
 *   column per locus (of the .map file)
 * @param locus_names the names of the loci (for error messages)
 * @param na_strings the alleles that are missing (e.g. "0")
 * @param n_threads the number of threads
 * @return a list with the first six fields of each line (`family_id`,
 *   `sample_id`, `paternal_id`, `maternal_id`, `sex` and `affection`, as
 *   strings), the first (`allele_ref`) and second (`allele_alt`) allele of
 *   each locus (NA if none), and the number of loci with no data
 *   (`n_no_data`), monomorphic (`n_monomorphic`) and multi-allelic
 *   (`n_multiallelic`, whose genotypes are all set to missing)
 *****************************************************************************/
// [[Rcpp::export]]
List ped_to_fbm(std::string filename, Environment BM,
                CharacterVector locus_names, CharacterVector na_strings,
                int n_threads = 1) {
  XPtr<FBM_RW> xp_bm = BM["address_rw"];
  unsigned char* fbm = static_cast<unsigned char*>(xp_bm->matrix());
  const size_t n_indiv = xp_bm->nrow();
  const size_t n_loci = xp_bm->ncol();
  if ((size_t) locus_names.size() != n_loci) {
    stop("'locus_names' should have one name per column of the FBM");
  }
  const unsigned char na_byte = 3;

  // the lines are read in chunks large enough to give each thread a few
  // lines to split (a line holds about 4 bytes per locus)
  size_t chunk_bytes = std::max((size_t) (1 << 24),
                                (size_t) n_threads * 2 * (4 * n_loci + 64));
  VcfLineReader reader(filename, n_threads, chunk_bytes, true);
  AlleleIds allele_ids(na_strings);

  // the state of each locus: the ids of its first and second allele (0 if
  // not found yet), and whether it is multi-allelic
  std::vector<uint16_t> allele_1(n_loci, 0);
  std::vector<uint16_t> allele_2(n_loci, 0);
  std::vector<unsigned char> multiallelic(n_loci, 0);

  // the first six fields of each individual
  StringBuffer fam[6];

  // the lines of the current chunk, the non blank ones (with their line
  // number in the file), what we find on each, and the ids of their alleles
  // (2 * n_loci per line, one line after the other)
  std::vector<char*> line_start;
  std::vector<size_t> line_len;
  std::vector<size_t> chunk_lines;
  std::vector<double> chunk_line_number;
  std::vector<PedLine> parsed;
  std::vector<uint16_t> ids;
  const size_t n_ids = 2 * n_loci;

  size_t n_read = 0; // the individuals read so far
  double line_number = 0; // the lines of the file read so far
  while (reader.next_chunk(line_start, line_len)) {
    checkUserInterrupt();
    // the lines with an individual (blank lines are skipped)
    chunk_lines.clear();
    chunk_line_number.clear();
    for (size_t l = 0; l < line_start.size(); l++) {
      line_number++;
      if (!is_blank_line(line_start[l], line_len[l])) {
        chunk_lines.push_back(l);
        chunk_line_number.push_back(line_number);
      }
    }
    const size_t n_lines = chunk_lines.size();
    if (n_lines == 0) {
      continue;
    }
    if (n_read + n_lines > n_indiv) {
      stop("the ped file has more individuals than when they were counted "
             "(did it change while being read?)");
    }
    parsed.resize(n_lines);
    if (ids.size() < n_lines * n_ids) {
      ids.resize(n_lines * n_ids);
    }

    // split the lines in parallel, turning each allele into its id
#pragma omp parallel for num_threads(n_threads) schedule(dynamic, 1)
    for (long b = 0; b < (long) n_lines; b++) {
      size_t l = chunk_lines[b];
      split_ped_line(line_start[l], line_len[l], n_loci, allele_ids,
                     &ids[b * n_ids], parsed[b]);
    }
    // check the lines (in order, so that the first error is reported), and
    // keep their first six fields
    for (size_t b = 0; b < n_lines; b++) {
      const PedLine &line = parsed[b];
      if (line.n_fields != 6 + n_ids) {
        char buf[160];
        std::snprintf(buf, sizeof(buf),
                      " has %.0f fields, but %.0f were expected (6, plus "
                      "two alleles for each of the %.0f loci of the map "
                      "file)", (double) line.n_fields,
                      (double) (6 + n_ids), (double) n_loci);
        stop("the line of " + individual_label(line, chunk_line_number[b]) +
          buf);
      }
      if (line.bad_locus >= 0) {
        stop("too many distinct alleles (longer than one character): "
               "the alleles of " +
               individual_label(line, chunk_line_number[b]) +
               " at locus " +
               std::string(locus_names[line.bad_locus]) +
               " could not be read");
      }
      for (int k = 0; k < 6; k++) {
        fam[k].push_back(line.fam_start[k], line.fam_len[k]);
      }
    }

    // code the genotypes, in parallel over blocks of loci (each locus is
    // coded from the individuals in order, from its state after the previous
    // chunks), writing them straight into the FBM (the genotypes of a locus
    // in this chunk are consecutive bytes of the backing file)
    const long block_loci = 256;
    const long n_blocks = (n_loci + block_loci - 1) / block_loci;
#pragma omp parallel for num_threads(n_threads) schedule(dynamic, 1)
    for (long blk = 0; blk < n_blocks; blk++) {
      size_t j_end = std::min(n_loci, (size_t) (blk + 1) * block_loci);
      for (size_t j = blk * block_loci; j < j_end; j++) {
        uint16_t a1 = allele_1[j];
        uint16_t a2 = allele_2[j];
        bool multi = multiallelic[j];
        unsigned char* out = fbm + j * n_indiv + n_read;
        const uint16_t* gt_ids = &ids[2 * j];
        for (size_t b = 0; b < n_lines; b++, gt_ids += n_ids) {
          bool one = false;
          bool two = false;
          for (int k = 0; k < 2; k++) {
            uint16_t ak = gt_ids[k];
            // a missing allele changes nothing
            if (ak == 0) {
              continue;
            }
            // the first allele found is the first allele of the locus
            if (a1 == 0) {
              a1 = ak;
            }
            if (ak == a1) {
              one = true;
            } else {
              // the first allele that differs from the first is the second
              if (a2 == 0) {
                a2 = ak;
              }
              if (ak == a2) {
                two = true;
              } else {
                // a third allele (in either position of the genotype)
                // makes the locus multi-allelic
                multi = true;
              }
            }
          }
          out[b] = one ? (two ? 1 : 0) : (two ? 2 : na_byte);
        }
        allele_1[j] = a1;
        allele_2[j] = a2;
        multiallelic[j] = multi;
      }
    }
    n_read += n_lines;
  }
  if (n_read != n_indiv) {
    stop("the ped file has fewer individuals than when they were counted "
           "(did it change while being read?)");
  }

  // the genotypes of multi-allelic loci are all missing; count loci with no
  // data, and monomorphic ones
  int n_no_data = 0;
  int n_monomorphic = 0;
  int n_multiallelic = 0;
  for (size_t j = 0; j < n_loci; j++) {
    if (multiallelic[j]) {
      std::memset(fbm + j * n_indiv, na_byte, n_indiv);
      n_multiallelic++;
    }
    if (allele_1[j] == 0) {
      n_no_data++;
    } else if (allele_2[j] == 0) {
      n_monomorphic++;
    }
  }

  // the alleles of the loci (as R strings)
  CharacterVector allele_ref(n_loci);
  CharacterVector allele_alt(n_loci);
  for (size_t j = 0; j < n_loci; j++) {
    allele_ref[j] = allele_ids.name(allele_1[j]);
    allele_alt[j] = allele_ids.name(allele_2[j]);
  }

  return List::create(
    _["family_id"] = fam[0].to_r(),
    _["sample_id"] = fam[1].to_r(),
    _["paternal_id"] = fam[2].to_r(),
    _["maternal_id"] = fam[3].to_r(),
    _["sex"] = fam[4].to_r(),
    _["affection"] = fam[5].to_r(),
    _["allele_ref"] = allele_ref,
    _["allele_alt"] = allele_alt,
    _["n_no_data"] = n_no_data,
    _["n_monomorphic"] = n_monomorphic,
    _["n_multiallelic"] = n_multiallelic
  );
}
