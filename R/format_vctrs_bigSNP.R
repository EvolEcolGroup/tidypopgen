#' @export
format.vctrs_bigSNP <- function(x, ..., width = 2) {
  # nothing to format if there are no individuals
  if (length(x) == 0) {
    return(character(0))
  }
  # prevent width from being greater than the number of loci
  if (width > nrow(show_loci(x))) {
    width <- nrow(show_loci(x))
  }
  # if there are no loci, there are no genotypes to show
  if (width == 0) {
    return(rep("[]", length(x)))
  }
  geno <- show_genotypes(x, loci_indices = seq_len(width))
  geno[is.na(geno)] <- "."
  geno <- paste0(
    "[",
    apply(geno, 1, function(row) {
      paste(row, collapse = ",")
    }), ",...]"
  )
  return(geno)
}
