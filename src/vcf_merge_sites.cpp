#include <Rcpp.h>
#include <vector>
#include <cstdint>

using namespace Rcpp;

// The sites of several merged VCFs (see vcf_merge_sites_cpp() below), with a
// hash table to find the sites at a given chromosome and position. The table
// (open addressing, with linear probing) holds the first site at each
// position; further sites at the same position (from loci that could not be
// merged into a biallelic SNP, or that are at the same position in a file)
// follow it in a chain (next_at_pos). Nothing is allocated per site but the
// elements of a few vectors, so that millions of sites are cheap.
class MergedSites {
public:
  // the chromosome (code), position and alleles (a character, 0 if missing)
  // of each site
  std::vector<int> chromosome;
  std::vector<int> position;
  std::vector<char> allele_ref;
  std::vector<char> allele_alt;
  // the file (0-based) and the locus of that file (0-based) where each site
  // was first found
  std::vector<int> first_file;
  std::vector<int> first_locus;
  // the number of files that have each site
  std::vector<int> n_files;
  // the last file that has each site (a site gets at most one locus of each
  // file: loci of the same file are never merged, as in bcftools merge)
  std::vector<int> last_file;
  // 0 if the site is kept; 1 if it is dropped because its alleles can not
  // form a biallelic SNP (bcftools would make it multiallelic); 2 if it is
  // dropped because its REF alleles differ between files
  std::vector<int> status;
  // the next site at the same position (-1 if none)
  std::vector<int> next_at_pos;

  MergedSites() : slots(1 << 10, -1) {}

  size_t size() const { return chromosome.size(); }

  // the first site at a chromosome and position (-1 if none)
  int find(int chr, int pos) const {
    size_t mask = slots.size() - 1;
    size_t i = hash(chr, pos) & mask;
    while (slots[i] >= 0) {
      int s = slots[i];
      if ((chromosome[s] == chr) && (position[s] == pos)) {
        return s;
      }
      i = (i + 1) & mask;
    }
    return -1;
  }

  // add a new site, found at locus `locus` of file `file`, after the last
  // site at its position (head is the first site at that position, -1 if
  // none); returns its index
  int add(int chr, int pos, char ref, char alt, int file, int locus,
          int head) {
    int s = size();
    chromosome.push_back(chr);
    position.push_back(pos);
    allele_ref.push_back(ref);
    allele_alt.push_back(alt);
    first_file.push_back(file);
    first_locus.push_back(locus);
    n_files.push_back(1);
    last_file.push_back(file);
    status.push_back(0);
    next_at_pos.push_back(-1);
    if (head >= 0) {
      // append the site to the chain of its position (in order of creation)
      int tail = head;
      while (next_at_pos[tail] >= 0) {
        tail = next_at_pos[tail];
      }
      next_at_pos[tail] = s;
    } else {
      // a new position: put it in the hash table
      size_t mask = slots.size() - 1;
      size_t i = hash(chr, pos) & mask;
      while (slots[i] >= 0) {
        i = (i + 1) & mask;
      }
      slots[i] = s;
      n_positions++;
      // keep the table less than 70% full (so that probes stay short)
      if (10 * n_positions > 7 * slots.size()) {
        grow();
      }
    }
    return s;
  }

private:
  std::vector<int> slots; // the first site at each position (-1 if empty)
  size_t n_positions = 0; // the number of positions in the table
  // hash of a chromosome and a position (the final mix of MurmurHash3, so
  // that all the bits of the slot depend on both)
  static uint32_t hash(int chr, int pos) {
    uint64_t h = ((uint64_t) (uint32_t) chr << 32) | (uint32_t) pos;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return (uint32_t) h;
  }
  // double the size of the hash table, and put back the first site of each
  // position
  void grow() {
    std::vector<int> old_slots(slots.size() * 2, -1);
    old_slots.swap(slots);
    size_t mask = slots.size() - 1;
    for (size_t j = 0; j < old_slots.size(); j++) {
      int s = old_slots[j];
      if (s < 0) {
        continue;
      }
      size_t i = hash(chromosome[s], position[s]) & mask;
      while (slots[i] >= 0) {
        i = (i + 1) & mask;
      }
      slots[i] = s;
    }
  }
};

// the allele of a character vector (of single characters) as a character, 0
// if it is missing (NA)
inline char allele_code(SEXP x, R_xlen_t i) {
  SEXP a = STRING_ELT(x, i);
  return (a == NA_STRING) ? 0 : CHAR(a)[0];
}

// whether two alleles agree: they are the same, or either is missing (an
// unobserved allele, e.g. the ALT "." of a monomorphic locus, can be any)
inline bool alleles_agree(char a, char b) {
  return (a == b) || (a == 0) || (b == 0);
}

// a character vector from alleles (NA for missing)
CharacterVector alleles_to_r(const std::vector<char> &x) {
  CharacterVector res(x.size());
  std::vector<SEXP> strings(256, R_NilValue);
  for (size_t i = 0; i < x.size(); i++) {
    unsigned char c = x[i];
    if (c == 0) {
      SET_STRING_ELT(res, i, NA_STRING);
      continue;
    }
    if (strings[c] == R_NilValue) {
      strings[c] = Rf_mkCharLen(&x[i], 1);
    }
    SET_STRING_ELT(res, i, strings[c]);
  }
  return res;
}

/******************************************************************************
 * Match the biallelic SNPs of several VCFs into merged sites
 *
 * The loci of the files are taken in order (all the loci of the first file,
 * then of the second, ...), and each locus is merged with a site at the same
 * chromosome and position, as in bcftools merge:
 * - a locus whose alleles agree with those of a site (the same REF and ALT,
 *   or a missing allele, which can be any) is merged with it (a missing
 *   allele of the site is then filled in from the locus);
 * - with ref_rule 2 ("swap"), a locus whose alleles agree with those of a
 *   site once swapped (its REF is the ALT of the site, and vice versa) is
 *   merged with it, and its genotypes will be swapped;
 * - a locus with a different REF from that of the site (both known) is a REF
 *   mismatch, which bcftools does not allow: with ref_rule 0 ("error") it is
 *   recorded (and gt_merge_vcfs() stops), with ref_rule 1 ("drop") the site
 *   is dropped (and the locus merged with it, so that later loci at the same
 *   position are dropped as well); with ref_rule 2 it is treated as below;
 * - other loci (the same REF, a different ALT: the merged site would be
 *   multiallelic, which a gen_tibble can not store) drop the site (as above),
 *   or, with split_alt (bcftools merge -m none), become a new site at the
 *   same position.
 * A locus is never merged with a site that already has a locus of the same
 * file (so loci at the same position in a file stay apart, as in bcftools),
 * and a locus with no site at its position becomes a new site.
 *
 * @param chromosome a list with, for each file, the chromosome (an integer
 *   code, the same in all files) of each biallelic SNP
 * @param position a list with, for each file, the position of each SNP
 * @param allele_ref a list with, for each file, the REF allele of each SNP
 *   (a character vector of single characters, NA if missing)
 * @param allele_alt as allele_ref, for the ALT allele
 * @param split_alt whether loci with different alleles become separate sites
 *   (rather than dropping their site)
 * @param ref_rule what to do with REF mismatches: 0 (record them, to stop
 *   with an error), 1 (drop the site), 2 (swap the alleles if they match)
 * @return a list with `site` (for each file, the 1-based site of each locus),
 *   `swap` (for each file, whether the alleles of each locus are swapped),
 *   and the vectors of each site (`chromosome`, `position`, `allele_ref`,
 *   `allele_alt`, `first_file` and `first_locus` (1-based), `n_files`,
 *   `status`), and `ref_mismatch`, the (1-based) site, file and locus of the
 *   first REF mismatches (with ref_rule 0), and their number
 */
// [[Rcpp::export]]
List vcf_merge_sites_cpp(const List& chromosome, const List& position,
                         const List& allele_ref, const List& allele_alt,
                         bool split_alt, int ref_rule) {
  int n_files = chromosome.size();
  if ((position.size() != n_files) || (allele_ref.size() != n_files) ||
      (allele_alt.size() != n_files)) {
    stop("chromosome, position, allele_ref and allele_alt should have one "
         "element per file");
  }
  if ((ref_rule < 0) || (ref_rule > 2)) {
    stop("ref_rule should be 0, 1 or 2");
  }
  MergedSites sites;
  List site_of_locus(n_files);
  List swap_of_locus(n_files);
  // the first REF mismatches (site, file and locus), and their number
  const int max_mismatches = 5;
  std::vector<int> mismatch_site, mismatch_file, mismatch_locus;
  double n_mismatches = 0;

  for (int f = 0; f < n_files; f++) {
    IntegerVector chr = chromosome[f];
    IntegerVector pos = position[f];
    CharacterVector ref = allele_ref[f];
    CharacterVector alt = allele_alt[f];
    R_xlen_t n_loci = chr.size();
    if ((pos.size() != n_loci) || (ref.size() != n_loci) ||
        (alt.size() != n_loci)) {
      stop("the loci of each file should have a chromosome, a position and "
           "two alleles");
    }
    IntegerVector site(n_loci);
    LogicalVector swap(n_loci);
    for (R_xlen_t k = 0; k < n_loci; k++) {
      int this_chr = chr[k];
      int this_pos = pos[k];
      char this_ref = allele_code(ref, k);
      char this_alt = allele_code(alt, k);
      int head = sites.find(this_chr, this_pos);
      // look for a site at this position that agrees with the locus (and
      // has no locus of this file yet), noting the first other site that
      // has no locus of this file (with which the locus conflicts)
      int match = -1;
      bool swapped = false;
      int conflict = -1;
      for (int s = head; s >= 0; s = sites.next_at_pos[s]) {
        if (sites.last_file[s] == f) {
          continue;
        }
        if (sites.status[s] == 0) {
          if (alleles_agree(this_ref, sites.allele_ref[s]) &&
              alleles_agree(this_alt, sites.allele_alt[s])) {
            match = s;
            break;
          }
          if ((ref_rule == 2) &&
              alleles_agree(this_ref, sites.allele_alt[s]) &&
              alleles_agree(this_alt, sites.allele_ref[s])) {
            match = s;
            swapped = true;
            break;
          }
        }
        if (conflict < 0) {
          conflict = s;
        }
      }
      int s;
      if (match >= 0) {
        // merge the locus with the site, filling in its missing alleles
        s = match;
        char &site_ref = sites.allele_ref[s];
        char &site_alt = sites.allele_alt[s];
        if (site_ref == 0) {
          site_ref = swapped ? this_alt : this_ref;
        }
        if (site_alt == 0) {
          site_alt = swapped ? this_ref : this_alt;
        }
      } else if (conflict >= 0) {
        s = conflict;
        // (a site that is already dropped takes the locus as it is)
        if (sites.status[s] == 0) {
          bool ref_differs = (this_ref != 0) && (sites.allele_ref[s] != 0) &&
            (this_ref != sites.allele_ref[s]);
          if (ref_differs && (ref_rule != 2)) {
            // a REF mismatch: drop the site (and record it, to stop with an
            // error if ref_rule is 0)
            sites.status[s] = 2;
            if (ref_rule == 0) {
              if (mismatch_site.size() < (size_t) max_mismatches) {
                mismatch_site.push_back(s + 1);
                mismatch_file.push_back(f + 1);
                mismatch_locus.push_back(k + 1);
              }
              n_mismatches++;
            }
          } else if (split_alt) {
            // alleles that can not be merged into a biallelic SNP: a new
            // site at the same position
            s = sites.add(this_chr, this_pos, this_ref, this_alt, f, k, head);
          } else {
            // a multiallelic site: drop it
            sites.status[s] = 1;
          }
        }
      } else {
        // a new site
        s = sites.add(this_chr, this_pos, this_ref, this_alt, f, k, head);
      }
      // (a new site already counts this file)
      if (sites.last_file[s] != f) {
        sites.last_file[s] = f;
        sites.n_files[s]++;
      }
      site[k] = s + 1;
      swap[k] = swapped;
    }
    site_of_locus[f] = site;
    swap_of_locus[f] = swap;
    Rcpp::checkUserInterrupt();
  }

  // the sites (with 1-based files and loci)
  IntegerVector first_file(sites.first_file.begin(), sites.first_file.end());
  IntegerVector first_locus(sites.first_locus.begin(),
                            sites.first_locus.end());
  first_file = first_file + 1;
  first_locus = first_locus + 1;
  return List::create(
    Named("site") = site_of_locus,
    Named("swap") = swap_of_locus,
    Named("chromosome") = IntegerVector(sites.chromosome.begin(),
                                        sites.chromosome.end()),
    Named("position") = IntegerVector(sites.position.begin(),
                                      sites.position.end()),
    Named("allele_ref") = alleles_to_r(sites.allele_ref),
    Named("allele_alt") = alleles_to_r(sites.allele_alt),
    Named("first_file") = first_file,
    Named("first_locus") = first_locus,
    Named("n_files") = IntegerVector(sites.n_files.begin(),
                                     sites.n_files.end()),
    Named("status") = IntegerVector(sites.status.begin(), sites.status.end()),
    Named("ref_mismatch") = List::create(
      Named("site") = mismatch_site,
      Named("file") = mismatch_file,
      Named("locus") = mismatch_locus,
      Named("n") = n_mismatches
    )
  );
}
