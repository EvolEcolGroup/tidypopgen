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
#' `position`, `ploidy`, `ind.names`. The `pop` slot is optional; if absent, the
#' returned gen_tibble will omit the population column.
#' @export
#' @examplesIf rlang::is_installed("adegenet")
#'
#' # Create a simple genlight object
#' x <- new("genlight",
#'   list(
#'     indiv1 = c(1, 1, 0, 1, 1, 0),
#'     indiv2 = c(2, 1, 1, 0, 0, 0)
#'   ),
#'   ploidy = c(2, 2),
#'   loc.names = paste0("locus", 1:6),
#'   chromosome = c("chr1", "chr1", "chr2", "chr2", "chr3", "chr3"),
#'   position = c(100, 200, 150, 250, 300, 400),
#'   loc.all = c("A/T", "C/G", "G/C", "A/T", "T/C", "G/A"),
#'   pop = c("pop1", "pop2")
#' )
#'
#'
#' file <- paste0(tempfile(), "gt_from_genlight")
#' # Convert to gen_tibble
#' new_gt <- gt_from_genlight(x, backingfile = file)
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
  indiv_meta <- data.frame(id = x@ind.names, stringsAsFactors = FALSE)
  if (!is.null(x@pop) && length(x@pop) == length(x@ind.names)) {
    indiv_meta$population <- x@pop
  }
  # extract genotypes as a matrix
  genotypes <- as.matrix(x)
  rownames(genotypes) <- NULL

  # create loci data.frame
  alleles <- x@loc.all
  # remove allele before / in every entry of alleles
  # adegenet checks that loc.all is in the form "A/T" or "G/C" etc
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
    ...,
    quiet = TRUE
  )
  return(gt)
}
