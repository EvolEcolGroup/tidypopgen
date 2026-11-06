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
  # recode the dosages only (not the population labels): hierfstat codes the
  # genotypes 0, 1 and 2 as 11, 12 and 22. The genotypes are decoded in
  # blocks of loci, and each locus is recoded directly into a column of the
  # data frame (recoding the whole matrix and then splitting it into the
  # columns of a data frame took more than five times the memory of the
  # result)
  codes <- c(11, 12, 22)
  n_loci <- count_loci(x)
  geno_cols <- vector("list", n_loci)
  # blocks of about 1 GB of decoded genotypes (smaller blocks use less
  # memory, but take longer)
  block_size <- bigstatsr::block_size(nrow(x))
  for (i_block in seq_len(ceiling(n_loci / block_size))) {
    loci_block <- seq(
      (i_block - 1) * block_size + 1,
      min(n_loci, i_block * block_size)
    )
    geno <- show_genotypes(x, loci_indices = loci_block)
    for (k in seq_along(loci_block)) {
      geno_cols[[loci_block[k]]] <- codes[geno[, k] + 1]
    }
  }
  names(geno_cols) <- loci_names(x)
  # the same data frame as data.frame(pop = ..., <matrix of genotypes>,
  # check.names = FALSE), built directly from the columns
  res <- c(list(pop = as.factor(x$population)), geno_cols)
  attr(res, "row.names") <- .set_row_names(nrow(x))
  class(res) <- "data.frame"
  res
}
