#include <Rcpp.h>
#include <vector>
#include <string>
#include <cstdint>

using namespace Rcpp;

// A site of several merged VCFs: its chromosome (code), position and alleles
// (a character, 0 if missing), the file (0-based) and the locus of that file
// (0-based) where it was first found, the number of files that have it, the
// last file that has it (a site gets at most one locus of each file: loci of
// the same file are never merged, as in bcftools merge), the next site at the
// same position (-1 if none), and its status: 0 if the site is kept; 1 if it
// is dropped because its alleles can not form a biallelic SNP (bcftools would
// make it multiallelic); 2 if it is dropped because its REF alleles differ
// between files. The fields of a site are stored together (32 bytes), so
// that looking up a site touches a single cache line.
struct MergedSite {
  int chromosome;
  int position;
  int first_file;
  int first_locus;
  int n_files;
  int last_file;
  int next_at_pos;
  char allele_ref;
  char allele_alt;
  char status;
};

// The sites of several merged VCFs, to which the loci of the files are added
// one file at a time (see merge_sites_new() below), with a hash table to find
// the sites at a given chromosome and position. The table (open addressing,
// with linear probing) holds the first site at each position; further sites
// at the same position (from loci that could not be merged into a biallelic
// SNP, or that are at the same position in a file) follow it in a chain
// (next_at_pos). For each file, it keeps the site of each locus, and whether
// its alleles are swapped (5 bytes per locus), until they are used to write
// the genotypes of the file (see merge_sites_file_cols()).
class MergedSites {
public:
  std::vector<MergedSite> sites;
  // the site (0-based) of each locus of each file, and whether its alleles
  // are swapped
  std::vector<std::vector<int> > site_of_locus;
  std::vector<std::vector<char> > swap_of_locus;
  std::vector<bool> released; // whether the loci of a file were released
  // the first REF mismatches: site, file and locus (0-based), and the alleles
  // of the site and of the locus; and the number of mismatches
  static const size_t max_mismatches = 5;
  std::vector<int> mismatch_site, mismatch_file, mismatch_locus;
  std::vector<char> mismatch_site_ref, mismatch_site_alt;
  std::vector<char> mismatch_locus_ref, mismatch_locus_alt;
  double n_mismatches;
  double n_swapped;

  MergedSites(bool split_alt, int ref_rule,
              const CharacterVector& missing_alleles)
    : n_mismatches(0), n_swapped(0), split_alt(split_alt),
      ref_rule(ref_rule), slots(1 << 10, -1), n_positions(0) {
    // missing alleles of a single character (alleles of biallelic SNPs are a
    // single character)
    for (int i = 0; i < 256; i++) {
      missing[i] = false;
    }
    for (R_xlen_t i = 0; i < missing_alleles.size(); i++) {
      SEXP a = STRING_ELT(missing_alleles, i);
      if ((a != NA_STRING) && (LENGTH(a) == 1)) {
        missing[(unsigned char) CHAR(a)[0]] = true;
      }
    }
  }

  // add the loci of the next file, matching them with the sites so far;
  // returns the (1-based) loci of the file that created new sites (in the
  // order of the sites)
  IntegerVector add_file(const IntegerVector& chr, const IntegerVector& pos,
                         const CharacterVector& ref,
                         const CharacterVector& alt) {
    R_xlen_t n_loci = chr.size();
    if ((pos.size() != n_loci) || (ref.size() != n_loci) ||
        (alt.size() != n_loci)) {
      stop("the loci of each file should have a chromosome, a position and "
           "two alleles");
    }
    int f = site_of_locus.size();
    std::vector<int> site(n_loci);
    std::vector<char> swap(n_loci, 0);
    std::vector<int> new_loci;
    const int* chr_ptr = chr.begin();
    const int* pos_ptr = pos.begin();
    for (R_xlen_t k = 0; k < n_loci; k++) {
      // with millions of sites, finding the site of a locus is dominated by
      // two cache misses (its slot in the hash table, then the site), so we
      // prefetch the slot of a locus a few loci ahead, and then its first
      // site (prefetching is only a hint: it does not change the results).
      // (This is written here rather than in a helper function, which GCC
      // drops as if it had no effect.)
#if defined(__GNUC__) || defined(__clang__)
      size_t mask = slots.size() - 1;
      if (k + 2 * prefetch_distance < n_loci) {
        R_xlen_t ahead = k + 2 * prefetch_distance;
        __builtin_prefetch(slots.data() +
                           (hash(chr_ptr[ahead], pos_ptr[ahead]) & mask));
      }
      if (k + prefetch_distance < n_loci) {
        R_xlen_t ahead = k + prefetch_distance;
        int s_ahead = slots[hash(chr_ptr[ahead], pos_ptr[ahead]) & mask];
        if (s_ahead >= 0) {
          __builtin_prefetch(sites.data() + s_ahead);
        }
      }
#endif
      int this_chr = chr_ptr[k];
      int this_pos = pos_ptr[k];
      char this_ref = allele_code(ref, k);
      char this_alt = allele_code(alt, k);
      int head = find(this_chr, this_pos);
      // look for a site at this position that agrees with the locus (and
      // has no locus of this file yet), noting the first other site that
      // has no locus of this file (with which the locus conflicts)
      int match = -1;
      bool swapped = false;
      int conflict = -1;
      for (int s = head; s >= 0; s = sites[s].next_at_pos) {
        const MergedSite &site_s = sites[s];
        if (site_s.last_file == f) {
          continue;
        }
        if (site_s.status == 0) {
          if (alleles_agree(this_ref, site_s.allele_ref) &&
              alleles_agree(this_alt, site_s.allele_alt)) {
            match = s;
            break;
          }
          if ((ref_rule == 2) &&
              alleles_agree(this_ref, site_s.allele_alt) &&
              alleles_agree(this_alt, site_s.allele_ref)) {
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
      bool new_site = false;
      if (match >= 0) {
        // merge the locus with the site, filling in its missing alleles
        s = match;
        MergedSite &site_s = sites[s];
        if (site_s.allele_ref == 0) {
          site_s.allele_ref = swapped ? this_alt : this_ref;
        }
        if (site_s.allele_alt == 0) {
          site_s.allele_alt = swapped ? this_ref : this_alt;
        }
      } else if (conflict >= 0) {
        s = conflict;
        // (a site that is already dropped takes the locus as it is)
        if (sites[s].status == 0) {
          bool ref_differs = (this_ref != 0) && (sites[s].allele_ref != 0) &&
            (this_ref != sites[s].allele_ref);
          if (ref_differs && (ref_rule != 2)) {
            // a REF mismatch: drop the site (and record it, to stop with an
            // error if ref_rule is 0)
            sites[s].status = 2;
            if (ref_rule == 0) {
              if (mismatch_site.size() < max_mismatches) {
                mismatch_site.push_back(s);
                mismatch_file.push_back(f);
                mismatch_locus.push_back(k);
                mismatch_site_ref.push_back(sites[s].allele_ref);
                mismatch_site_alt.push_back(sites[s].allele_alt);
                mismatch_locus_ref.push_back(this_ref);
                mismatch_locus_alt.push_back(this_alt);
              }
              n_mismatches++;
            }
          } else if (split_alt) {
            // alleles that can not be merged into a biallelic SNP: a new
            // site at the same position
            s = add(this_chr, this_pos, this_ref, this_alt, f, k, head);
            new_site = true;
          } else {
            // a multiallelic site: drop it
            sites[s].status = 1;
          }
        }
      } else {
        // a new site
        s = add(this_chr, this_pos, this_ref, this_alt, f, k, head);
        new_site = true;
      }
      if (new_site) {
        new_loci.push_back(k + 1);
      }
      // (a new site already counts this file)
      if (sites[s].last_file != f) {
        sites[s].last_file = f;
        sites[s].n_files++;
      }
      site[k] = s;
      if (swapped) {
        swap[k] = 1;
        n_swapped++;
      }
    }
    site_of_locus.push_back(std::vector<int>());
    site_of_locus.back().swap(site);
    swap_of_locus.push_back(std::vector<char>());
    swap_of_locus.back().swap(swap);
    released.push_back(false);
    return IntegerVector(new_loci.begin(), new_loci.end());
  }

private:
  bool split_alt;
  int ref_rule;
  bool missing[256];
  std::vector<int> slots; // the first site at each position (-1 if empty)
  size_t n_positions; // the number of positions in the table

  // the allele of a character vector (of single characters) as a character,
  // 0 if it is missing (NA or a missing allele)
  char allele_code(const CharacterVector& x, R_xlen_t i) const {
    SEXP a = STRING_ELT(x, i);
    if (a == NA_STRING) {
      return 0;
    }
    unsigned char c = CHAR(a)[0];
    return missing[c] ? 0 : (char) c;
  }

  // whether two alleles agree: they are the same, or either is missing (an
  // unobserved allele, e.g. the ALT "." of a monomorphic locus, can be any)
  static bool alleles_agree(char a, char b) {
    return (a == b) || (a == 0) || (b == 0);
  }

  // how many loci ahead we prefetch the first site at a position (and twice
  // as many for its slot in the hash table)
  static const R_xlen_t prefetch_distance = 8;


  // the first site at a chromosome and position (-1 if none)
  int find(int chr, int pos) const {
    size_t mask = slots.size() - 1;
    size_t i = hash(chr, pos) & mask;
    while (slots[i] >= 0) {
      const MergedSite &site_s = sites[slots[i]];
      if ((site_s.chromosome == chr) && (site_s.position == pos)) {
        return slots[i];
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
    int s = sites.size();
    MergedSite site_s = {chr, pos, file, locus, 1, file, -1, ref, alt, 0};
    sites.push_back(site_s);
    if (head >= 0) {
      // append the site to the chain of its position (in order of creation)
      int tail = head;
      while (sites[tail].next_at_pos >= 0) {
        tail = sites[tail].next_at_pos;
      }
      sites[tail].next_at_pos = s;
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
      size_t i = hash(sites[s].chromosome, sites[s].position) & mask;
      while (slots[i] >= 0) {
        i = (i + 1) & mask;
      }
      slots[i] = s;
    }
  }
};

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

// the sites of an external pointer created by merge_sites_new()
MergedSites* get_merged_sites(SEXP sites) {
  XPtr<MergedSites> ptr(sites);
  if (ptr.get() == NULL) {
    stop("the merged sites are no longer valid");
  }
  return ptr.get();
}

/******************************************************************************
 * Match the biallelic SNPs of several VCFs into merged sites
 *
 * The loci of the files are added one file at a time (with
 * merge_sites_add_file(), all the loci of the first file, then of the
 * second, ...), and each locus is merged with a site at the same chromosome
 * and position, as in bcftools merge:
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
 * Only the sites, and the site of each locus (with whether its alleles are
 * swapped), are kept: the loci of each file can be dropped once they have
 * been added, so that the memory used does not grow with the loci of all the
 * files at once.
 *
 * @param split_alt whether loci with different alleles become separate sites
 *   (rather than dropping their site)
 * @param ref_rule what to do with REF mismatches: 0 (record them, to stop
 *   with an error), 1 (drop the site), 2 (swap the alleles if they match)
 * @param missing_alleles the alleles that are missing (besides NA)
 * @return an external pointer to the merged sites
 */
// [[Rcpp::export]]
SEXP merge_sites_new(bool split_alt, int ref_rule,
                     const CharacterVector& missing_alleles) {
  if ((ref_rule < 0) || (ref_rule > 2)) {
    stop("ref_rule should be 0, 1 or 2");
  }
  XPtr<MergedSites> ptr(new MergedSites(split_alt, ref_rule, missing_alleles),
                        true);
  return ptr;
}

/*
 * Add the biallelic SNPs of the next file to the merged sites
 *
 * @param sites the merged sites (from merge_sites_new())
 * @param chromosome the chromosome (an integer code, the same in all files)
 *   of each biallelic SNP of the file
 * @param position the position of each SNP
 * @param allele_ref the REF allele of each SNP (a character vector of single
 *   characters)
 * @param allele_alt as allele_ref, for the ALT allele
 * @return the (1-based) SNPs of the file that created new sites, in the order
 *   of the sites (so that the sites can take the names of these SNPs)
 */
// [[Rcpp::export]]
IntegerVector merge_sites_add_file(SEXP sites, const IntegerVector& chromosome,
                                   const IntegerVector& position,
                                   const CharacterVector& allele_ref,
                                   const CharacterVector& allele_alt) {
  return get_merged_sites(sites)->add_file(chromosome, position, allele_ref,
                                           allele_alt);
}

/*
 * The merged sites
 *
 * @param sites the merged sites (from merge_sites_new())
 * @return a list with the vectors of each site (`chromosome`, `position`,
 *   `allele_ref`, `allele_alt`, `first_file` and `first_locus` (1-based),
 *   `n_files`, `status`), `n_swapped`, the number of loci whose alleles are
 *   swapped, and `ref_mismatch`, the (1-based) site, file and locus of the
 *   first REF mismatches (with ref_rule 0), with the alleles of the site and
 *   of the locus, and their number
 */
// [[Rcpp::export]]
List merge_sites_result(SEXP sites) {
  MergedSites* ms = get_merged_sites(sites);
  size_t n_sites = ms->sites.size();
  IntegerVector chromosome(n_sites), position(n_sites);
  IntegerVector first_file(n_sites), first_locus(n_sites);
  IntegerVector n_files(n_sites), status(n_sites);
  std::vector<char> allele_ref(n_sites), allele_alt(n_sites);
  for (size_t s = 0; s < n_sites; s++) {
    const MergedSite &site_s = ms->sites[s];
    chromosome[s] = site_s.chromosome;
    position[s] = site_s.position;
    first_file[s] = site_s.first_file + 1;
    first_locus[s] = site_s.first_locus + 1;
    n_files[s] = site_s.n_files;
    status[s] = site_s.status;
    allele_ref[s] = site_s.allele_ref;
    allele_alt[s] = site_s.allele_alt;
  }
  IntegerVector mismatch_site(ms->mismatch_site.begin(),
                              ms->mismatch_site.end());
  IntegerVector mismatch_file(ms->mismatch_file.begin(),
                              ms->mismatch_file.end());
  IntegerVector mismatch_locus(ms->mismatch_locus.begin(),
                               ms->mismatch_locus.end());
  return List::create(
    Named("chromosome") = chromosome,
    Named("position") = position,
    Named("allele_ref") = alleles_to_r(allele_ref),
    Named("allele_alt") = alleles_to_r(allele_alt),
    Named("first_file") = first_file,
    Named("first_locus") = first_locus,
    Named("n_files") = n_files,
    Named("status") = status,
    Named("n_swapped") = ms->n_swapped,
    Named("ref_mismatch") = List::create(
      Named("site") = mismatch_site + 1,
      Named("file") = mismatch_file + 1,
      Named("locus") = mismatch_locus + 1,
      Named("site_ref") = alleles_to_r(ms->mismatch_site_ref),
      Named("site_alt") = alleles_to_r(ms->mismatch_site_alt),
      Named("locus_ref") = alleles_to_r(ms->mismatch_locus_ref),
      Named("locus_alt") = alleles_to_r(ms->mismatch_locus_alt),
      Named("n") = ms->n_mismatches
    )
  );
}

/*
 * The columns of the merged FBM of the biallelic SNPs of a file
 *
 * The sites of the SNPs of the file are then released (each file is used
 * once, to write its genotypes).
 *
 * @param sites the merged sites (from merge_sites_new())
 * @param file the (1-based) file
 * @param site_col the (0-based) column of each site in the merged FBM (-1 if
 *   the site is not kept)
 * @return a list with `col`, the column of each SNP of the file, and `swap`,
 *   whether its alleles are swapped (an empty vector if none is)
 */
// [[Rcpp::export]]
List merge_sites_file_cols(SEXP sites, int file,
                           const IntegerVector& site_col) {
  MergedSites* ms = get_merged_sites(sites);
  if ((file < 1) || ((size_t) file > ms->site_of_locus.size()) ||
      ms->released[file - 1]) {
    stop("the sites of this file are not available");
  }
  if ((size_t) site_col.size() != ms->sites.size()) {
    stop("'site_col' should have one value per site");
  }
  std::vector<int> &site = ms->site_of_locus[file - 1];
  std::vector<char> &swap = ms->swap_of_locus[file - 1];
  IntegerVector col(site.size());
  bool any_swap = false;
  for (size_t k = 0; k < site.size(); k++) {
    col[k] = site_col[site[k]];
    any_swap = any_swap || swap[k];
  }
  LogicalVector swap_r(any_swap ? swap.size() : 0);
  if (any_swap) {
    for (size_t k = 0; k < swap.size(); k++) {
      swap_r[k] = swap[k];
    }
  }
  // release the sites of the loci of this file
  std::vector<int>().swap(site);
  std::vector<char>().swap(swap);
  ms->released[file - 1] = true;
  return List::create(Named("col") = col, Named("swap") = swap_r);
}
