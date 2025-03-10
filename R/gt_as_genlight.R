#' Convert a `gen_tibble` to a `genlight` object from `adegenet`
#'
#' This function converts a `gen_tibble` to a `genlight` object from `adegenet`
#'
#' @param x a [`gen_tibble`], with population coded as 'population'
#' @returns a `genlight` object
#' @export
#'

gt_as_genlight <- function(x) {
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop(
      "to use this function, first install package 'adegenet' with\n",
      "install.packages('adegenet')"
    )
  }
  stopifnot_diploid(x) # currently we only support diploid data
  test_genlight <- methods::new("genlight", #nolint
    gen = show_genotypes(x),
    ploidy = 2, # TODO update this when ploidy is implemented
    ind.names = x$id,
    pop = x$population,
    loc.names = loci_names(x)
  )
}
