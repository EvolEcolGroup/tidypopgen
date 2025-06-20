#' Convert a `gen_tibble` to a `genlight` object from `adegenet`
#'
#' This function converts a `gen_tibble` to a `genlight` object from `adegenet`
#'
#' @param x a [`gen_tibble`], with population coded as 'population'
#' @returns a `genlight` object
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Convert to genind
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
  test_genlight <- methods::new( # nolint
    "genlight", # nolint
    gen = show_genotypes(x),
    ploidy = 2, # TODO update this when ploidy is implemented
    ind.names = x$id,
    pop = x$population,
    loc.names = loci_names(x)
  )
}
