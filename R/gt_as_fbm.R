#' Convert a `gen_tibble` to a file backed matrix `FBM` from `bigstatsr`
#'
#' This function converts a `gen_tibble` object to a file-backed matrix (`FBM`)
#' from the `bigstatsr` package. As `gen_tibble` objects use an `FBM`
#' internally, this simply returns the appropriate `FBM` object. Note that, for
#' tibble with imputed data, the genotype coding will depend on the impute
#' setting of the `gen_tibble` object (this can be checked with 
#' [gt_uses_imputed()] and changed with [gt_set_imputed()]).
#'
#' @param x a [`gen_tibble`]
#' @returns a [`FBM`] object from the `bigstatsr` package.
#' @export
#' @examples
#' # Create a gen_tibble
#' example_gt <- load_example_gt("gen_tbl")
#' # Convert the gen_tibble to a file-backed matrix
#' example_fbm <- gt_as_fbm(example_gt)
#' # Check the class of the resulting object
#' class(example_fbm)

gt_as_fbm <- function(x) {
  # Check if the input is a gen_tibble
  if (!inherits(x, "gen_tibble")) {
    stop("Input must be a gen_tibble object.")
  }
  
  # we use the internal function
  return(.gt_as_fbm(x))
}