#' Convert a `genlight` object from adegenet to a `gen_tibble`
#'
#' This function converts a `genlight` object from the `adegenet` package to a
#' `gen_tibble` object
#'
#' @param x A `genlight` object
#' @param backingfile the path, including the file name without extension, for
#'   backing files used to store the data (they will be given a .bk and .rds
#'   automatically). If `NULL` (default), backing files are placed in the
#'   temporary directory.
#' @param ... Additional arguments passed to gen_tibble().
#' @return A `gen_tibble` object
#' @details
#' - Currently supports diploid `genlight` objects only (all values in
#' `@ploidy` must be 2).
#' - Requires non-missing slots: `loc.names`, `n.loc`, `loc.all`, `chromosome`,
#' `position`, `ploidy`, `ind.names`, `pop`.
#' @export
#'
gt_from_genlight <- function(x, backingfile = NULL, ...) {
  if (is.null(backingfile)) {
    backingfile <- tempfile()
  }
  # check x contains the slots we need
  required_slots <- c(
    "loc.names", "n.loc", "loc.all", "chromosome", "position",
    "ploidy", "ind.names", "pop"
  )
  if (any(sapply(required_slots, function(slot) is.null(methods::slot(x, slot))))) { # nolint
    null_slots <-
      sapply(required_slots, function(slot) is.null(methods::slot(x, slot)))
    stop(
      "The genlight object has one or more required slots that are NULL: ",
      paste(required_slots[null_slots], collapse = ", ")
    )
  }
  # check @ploidy is 2
  if (any(x@ploidy != 2)) {
    stop(
      "Currently only diploid genlight objects are supported (i.e., all ",
      "values in the @ploidy slot must be 2)."
    )
  }

  # create indiv_meta data.frame
  indiv_meta <- data.frame(
    id = x@ind.names,
    population = x@pop
  )
  # extract genotypes as a matrix
  genotypes <- as.matrix(x)
  rownames(genotypes) <- NULL

  # create loci data.frame
  alleles <- x@loc.all
  if (any(!grepl("/", alleles, fixed = TRUE))) {
    stop("All loci must be biallelic and encoded like 'A/T'.")
  }
  # remove allele before / in every entry of alleles
  allele_ref <- gsub("/.*", "", alleles)
  allele_alt <- gsub(".*?/", "", alleles)
  loci <- data.frame(
    name = x@loc.names,
    chromosome = x@chromosome,
    position = x@position,
    genetic_dist = as.double(rep(0, x@n.loc)),
    allele_ref = allele_ref,
    allele_alt = allele_alt
  )

  gt <- gen_tibble(
    x = genotypes,
    loci = loci,
    indiv_meta = indiv_meta,
    backingfile = backingfile,
    quiet = TRUE,
    ...
  )
  return(gt)
}
