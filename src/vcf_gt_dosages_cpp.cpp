#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************/
// Dosages of the genotypes (GT strings, e.g. "0/1", "1|1", "./.", or
// "0/1/1/0" for polyploids) of a vcfR chunk, used by vcf_to_fbm_vcfR().
//
// This replaces apply(gt, 2, poly_indiv_dosage), which called strsplit() and
// an R function (poly_genotype_dosage()) for every genotype, with the same
// rules: alleles are separated by "/" or "|" (a final separator adds no
// allele, as with strsplit()); a genotype with more than max_ploidy alleles
// is an error; a genotype is missing (max_ploidy + 1) if it is NA, empty or
// has a missing allele ("."); otherwise its dosage is the sum of its allele
// indices, which must be smaller than max_ploidy + 1.
//
// gt: a character matrix of genotypes, with loci as rows and individuals as
//   columns (as returned by vcfR::extract.gt())
// Returns a list with `dosages`, a raw matrix with individuals as rows and
// loci as columns (as in the FBM), and `error`: 0 if all genotypes are
// valid, 1 if a genotype has more than max_ploidy alleles or a dosage larger
// than max_ploidy, and 2 if an allele is not a number (in which case the
// dosages are not complete).

// [[Rcpp::export]]
List vcf_gt_dosages_cpp(const CharacterMatrix& gt, int max_ploidy) {
  size_t n_loci = gt.nrow();
  size_t n_indiv = gt.ncol();
  RawMatrix dosages(n_indiv, n_loci);
  const unsigned char missing = (unsigned char) (max_ploidy + 1);
  int error = 0;
  for (size_t i = 0; i < n_indiv && error == 0; i++) {
    for (size_t j = 0; j < n_loci; j++) {
      SEXP g = STRING_ELT(gt, i * n_loci + j);
      if (g == NA_STRING || CHAR(g)[0] == '\0') {
        dosages(i, j) = missing;
        continue;
      }
      const char* s = CHAR(g);
      // read the alleles: count them, add up their indices, and flag
      // missing (".") and invalid alleles
      int n_alleles = 0;
      long dosage = 0;
      bool is_missing = false;
      bool is_invalid = false;
      const char* p = s;
      while (true) {
        // the allele from p to the next separator (or the end)
        const char* q = p;
        while (*q != '\0' && *q != '/' && *q != '|') q++;
        // a final separator does not add an empty allele (as strsplit())
        if (q == p && *q == '\0' && n_alleles > 0) break;
        n_alleles++;
        if (q - p == 1 && *p == '.') {
          is_missing = true;
        } else if (q == p) {
          // an empty allele (e.g. "0//1")
          is_invalid = true;
        } else {
          long value = 0;
          for (const char* c = p; c < q; c++) {
            if (*c < '0' || *c > '9') {
              is_invalid = true;
              break;
            }
            value = value * 10 + (*c - '0');
            if (value > 255) value = 256; // large enough to be an error
          }
          dosage += value;
        }
        if (*q == '\0') break;
        p = q + 1;
      }
      if (n_alleles > max_ploidy) {
        error = 1;
        break;
      }
      if (is_missing) {
        dosages(i, j) = missing;
      } else if (is_invalid) {
        error = 2;
        break;
      } else if (dosage >= max_ploidy + 1) {
        error = 1;
        break;
      } else {
        dosages(i, j) = (unsigned char) dosage;
      }
    }
  }
  return List::create(_["dosages"] = dosages, _["error"] = error);
}

/******************************************************************************/
