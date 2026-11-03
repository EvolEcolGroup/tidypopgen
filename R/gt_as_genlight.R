#' Convert a `gen_tibble` to a `genlight` object from `adegenet`
#'
#' This function converts a `gen_tibble` to a `genlight` object from `adegenet`
#'
#' The chromosome and position of each locus are stored in the `genlight`
#' object, and the alleles are stored as "ref/alt" in `loc.all` (`genlight`
#' counts the second allele, i.e. the alternate allele). As `genlight` only
#' supports single letter alleles, `loc.all` is not set (with a warning) if any
#' locus has a missing or multi-character allele.
#'
#' @param x a [`gen_tibble`], with population coded as 'population'
#' @returns a `genlight` object
#' @export
#' @examplesIf rlang::is_installed("adegenet")
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Convert to genlight
#' gt_genlight <- example_gt %>% gt_as_genlight()
#'
#' # Check object class
#' class(gt_genlight)
#'
gt_as_genlight <- function(x) {
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop(
      "to use this function, first install package 'adegenet' with\n",
      "install.packages('adegenet')"
    )
  }
  stopifnot_diploid(x) # currently we only support diploid data
  loci <- show_loci(x)
  # genlight counts the second allele in loc.all, i.e. the alternate allele
  loc_all <- paste0(loci$allele_ref, "/", loci$allele_alt)
  # adegenet only accepts single letter alleles (e.g. "a/t") for all loci
  if (!all(grepl("^[[:alpha:]]{1}/[[:alpha:]]{1}$", loc_all))) {
    warning(
      "some loci have missing or multi-character alleles, which genlight ",
      "does not support; the alleles (loc.all) were not stored"
    )
    loc_all <- NULL
  }
  test_genlight <- methods::new( # nolint
    "genlight", # nolint
    gen = show_genotypes(x),
    ploidy = 2, # TODO update this when ploidy is implemented
    ind.names = x$id,
    pop = x$population,
    loc.names = loci_names(x),
    loc.all = loc_all,
    chromosome = loci$chromosome,
    position = loci$position
  )
}
