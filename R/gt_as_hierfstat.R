#' Convert a `gen_tibble` to a data.frame compatible with `hierfstat`
#'
#' This function converts a `gen_tibble` to a data.frame formatted
#' to be used by `hierfstat` functions.
#'
#' @note `hierfstat` only supports haploid or diploid data, so this function
#'   only works on diploid (or pseudohaploid) `gen_tibble` objects.
#'   Pseudohaploid data is simply treated as regular diploid genotype counts
#'   (dosage 0/2), the same convention used by other tools (e.g. PLINK).
#' @param x a [`gen_tibble`], with population coded as 'population'
#' @returns a data.frame with a column 'pop' and further column representing
#' the genotypes (with alleles recoded as 1 and 2)
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Convert to hierfstat format
#' gt_hierfstat <- example_gt %>% gt_as_hierfstat()
#'
#' # Check object class
#' class(gt_hierfstat)
gt_as_hierfstat <- function(x) {
  stopifnot_dip_pseudo(x)
  # recode the dosages only (not the population labels)
  geno <- show_genotypes(x)
  geno[] <- c(11, 12, 22)[geno + 1]
  colnames(geno) <- loci_names(x)
  data.frame(pop = as.factor(x$population), geno, check.names = FALSE)
}
