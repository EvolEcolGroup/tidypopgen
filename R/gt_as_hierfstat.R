#' Convert a `gen_tibble` to a data.frame compatible with `hierfstat`
#'
#' This function converts a `gen_tibble` to a data.frame formatted
#' to be used by `hierfstat` functions.
#'
#' @param x a [`gen_tibble`], with population coded as 'population'
#' @returns a data.frame with a column 'pop' and further column representing
#' the genotypes (with alleles recoded as 1 and 2)
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Convert to genind
#' gt_hierfstat <- example_gt %>% gt_as_hierfstat()
#'
#' # Check object class
#' class(gt_hierfstat)
gt_as_hierfstat <- function(x) {
  hier_df <- data.frame(pop = as.factor(x$population), show_genotypes(x))
  hier_df[hier_df == 0] <- 11
  hier_df[hier_df == 1] <- 12
  hier_df[hier_df == 2] <- 22
  hier_df
}
