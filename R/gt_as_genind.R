#' Convert a `gen_tibble` to a `genind` object from `adegenet`
#'
#' This function converts a `gen_tibble` to a `genind` object from `adegenet`.
#' `genind` objects support arbitrary ploidy (including mixed-ploidy data
#' across individuals), so this works on diploid, polyploid and mixed-ploidy
#' `gen_tibble` objects alike.
#'
#' @note pseudohaploid individuals are treated as diploid (ploidy 2) using
#'   their stored genotype counts (dosage 0/2) directly, the same convention
#'   used by other tools (e.g. PLINK) when handling pseudohaploid data.
#' @param x a [`gen_tibble`], with population coded as 'population'
#' @returns a `genind` object
#' @export
#' @examplesIf rlang::is_installed("adegenet")
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Convert to genind
#' gt_genind <- example_gt %>% gt_as_genind()
#'
#' # Check object class
#' class(gt_genind)
gt_as_genind <- function(x) {
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop(
      "to use this function, first install package 'adegenet' with\n",
      "install.packages('adegenet')"
    )
  }
  geno_mat <- show_genotypes(x)
  # pseudohaploid individuals have their dosage stored on the doubled,
  # diploid-like scale (0/2, never 1); treat them as ploidy 2 when decoding
  # dosage into allele codes, matching every other diploid/pseudohaploid
  # dispatch in the package
  geno_ploidy <- indiv_ploidy(x)
  if (is_pseudohaploid(x)) {
    geno_ploidy[geno_ploidy == 1] <- 2
  }

  # build a genotype string per individual/locus (ncode = 1 digit per
  # allele, as required by adegenet::df2genind()): `dosage` copies of
  # allele "2" (alt) and `ploidy - dosage` copies of allele "1" (ref)
  dosage_to_genotype <- function(dosage, ploidy) {
    if (is.na(dosage)) {
      return(NA_character_)
    }
    paste0(strrep("1", ploidy - dosage), strrep("2", dosage))
  }
  df_for_genind <- matrix(
    mapply(
      dosage_to_genotype,
      dosage = as.vector(geno_mat),
      ploidy = rep(geno_ploidy, times = ncol(geno_mat))
    ),
    nrow = nrow(geno_mat),
    ncol = ncol(geno_mat)
  )

  adegenet::df2genind(
    X = df_for_genind,
    ind.names = x$id,
    pop = x$population,
    ploidy = geno_ploidy,
    ncode = 1,
    loc.names = loci_names(x)
  )
}
