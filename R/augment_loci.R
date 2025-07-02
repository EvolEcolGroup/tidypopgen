#' Augment the loci table with information from a analysis object
#'
#' `augment_loci` add columns to the loci table of a `gen_tibble` related to
#' information from a given analysis.
#' @param x  An object returned by one of the `gt_` functions (e.g. [gt_pca()]).
#' @param data the `gen_tibble` used to run the PCA.
#' @param ... Additional parameters passed to the individual methods.
#' @return A loci tibble with additional columns. If `data` is missing, a tibble
#'   of the information, with a column `.rownames` giving the loci names.
#' @export
#' @examples
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' lobsters <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Remove monomorphic loci and impute
#' lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
#' lobsters <- gt_impute_simple(lobsters, method = "mode")
#'
#' # Create PCA
#' pca <- gt_pca_partialSVD(lobsters)
#'
#' # Augment the gen_tibble with the PCA scores
#' augment_loci(pca, data = lobsters)
#'
augment_loci <- function(x, data, ...) {
  UseMethod("augment_loci", x)
}
