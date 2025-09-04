#' @export
format.vctrs_bigSNP <- function(x, ..., width = 2) {
  # prevent width from being greater than the number of loci
  if (width > nrow(show_loci(x))) {
    width <- nrow(show_loci(x))
  }
  geno <- show_genotypes(x, loci_indices = 1:width)
  geno[is.na(geno)] <- "."
  geno <- paste0(
    "[",
    apply(geno, 1, function(row) {
      paste(row, collapse = ",")
    }), ",...]"
  )
  return(geno)
}
