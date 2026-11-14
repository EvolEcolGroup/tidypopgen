#' Individual inbreeding coefficient
#'
#' This function calculates the inbreeding coefficient for each individual based
#' on the beta estimate from Weir and Goudet (2017).
#'
#' `indiv_inbreeding()` is an S3 generic, with methods for [`gen_tibble`]
#' objects (grouped or not) and for their `genotypes` column. Other packages
#' can add methods for their own classes (e.g. tibbles of genotype
#' likelihoods).
#'
#' @references Weir, BS and Goudet J (2017) A Unified Characterization of
#'   Population Structure and Relatedness. Genetics (2017) 206:2085
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`]. The method for the
#'   `genotypes` column needs `allele_sharing_mat`.
#' @param method currently only "WG17" (for Weir and Goudet 2017).
#' @param allele_sharing_mat optional and only relevant for "WG17", the matrix
#'   of Allele Sharing returned by [pairwise_allele_sharing()] with
#'   `as_matrix=TRUE`. As a number of statistics can be derived from the Allele
#'   Sharing matrix, it is sometimes more efficient to pre-compute this matrix.
#'   It is not possible to use this with grouped tibbles.
#' @param ... other arguments passed to specific methods, currently unused.
#' @return a numeric vector of inbreeding coefficients. For a grouped
#'   [`gen_tibble`], a list of such vectors, one per group (the coefficients
#'   are computed within each group).
#' @rdname indiv_inbreeding
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% indiv_inbreeding(method = "WG17")
#'
indiv_inbreeding <- function(.x, ...) {
  UseMethod("indiv_inbreeding", .x)
}

#' @export
#' @rdname indiv_inbreeding
indiv_inbreeding.gen_tbl <- function(
  .x,
  method = c("WG17"),
  allele_sharing_mat = NULL,
  ...
) {
  stopifnot_has_genotypes(.x)
  # the matrix of allele sharing is computed from the whole gen_tibble (it is
  # labelled with its ids), then the genotypes column is handed over to its
  # method
  if (is.null(allele_sharing_mat)) {
    allele_sharing_mat <- pairwise_allele_sharing(.x, as_matrix = TRUE)
  }
  indiv_inbreeding(
    .x$genotypes,
    method = method,
    allele_sharing_mat = allele_sharing_mat,
    ...
  )
}

#' @export
#' @rdname indiv_inbreeding
indiv_inbreeding.vctrs_bigSNP <- function(
  .x,
  method = c("WG17"),
  allele_sharing_mat = NULL,
  ...
) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  if (is.null(allele_sharing_mat)) {
    # TODO deal with the fact that allele sharing can only be computed
    # for the gen_tibble
    # TODO modify allele sharing code to be able to create a matrix just
    # from the genotype column
    stop(paste(
      "this method has to be applied to a gen_tibble, not just",
      "the genotype column"
    ))
  }

  # taken from beta.dosage in hierfstat
  mii <- diag(allele_sharing_mat) # diagonal of the matrix
  diag(allele_sharing_mat) <- NA
  mb <- mean(allele_sharing_mat, na.rm = TRUE)
  indiv_inb <- ((mii * 2 - 1) - mb) / (1 - mb)
  return(indiv_inb)
}

#' @export
#' @rdname indiv_inbreeding
indiv_inbreeding.grouped_gen_tbl <- function(
  .x,
  method = c("WG17"),
  allele_sharing_mat = NULL,
  ...
) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  if (!is.null(allele_sharing_mat)) {
    stop("allele_sharing_mat can not be provided for grouped_df objects")
  }
  group_map(
    .x,
    .f = ~ indiv_inbreeding(.x, method = method, allele_sharing_mat = NULL)
  )
}

# a grouped data frame that has lost the class of its tibble (e.g. a grouped
# gen_tibble returned by summarise()): restore the class and dispatch again
#' @export
#' @noRd
indiv_inbreeding.grouped_df <- function(.x, ...) {
  indiv_inbreeding(restore_genotypes_tbl(.x, "indiv_inbreeding"), ...)
}

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
indiv_inbreeding.data.frame <- function(.x, ...) {
  indiv_inbreeding(restore_genotypes_tbl(.x, "indiv_inbreeding"), ...)
}
