#' Individual inbreeding coefficient
#'
#' This function calculates the inbreeding coefficient for each individual based
#' on the beta estimate from Weir and Goudet (2017).
#'
#' @references Weir, BS and Goudet J (2017) A Unified Characterization of
#'   Population Structure and Relatedness. Genetics (2017) 206:2085
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param method currently only "WG17" (for Weir and Goudet 2017).
#' @param allele_sharing_mat optional and only relevant for "WG17", the matrix
#'   of Allele Sharing returned by [pairwise_allele_sharing()] with
#'   `as_matrix=TRUE`. As a number of statistics can be derived from the Allele
#'   Sharing matrix, it is sometimes more efficient to pre-compute this matrix.
#'   It is not possible to use this with grouped tibbles.
#' @param ... currently unused.
#' @return a numeric vector of inbreeding coefficients.
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% indiv_inbreeding(method = "WG17")
#'
indiv_inbreeding <- function(
    .x,
    method = c("WG17"),
    allele_sharing_mat = NULL,
    ...) {
  UseMethod("indiv_inbreeding", .x)
}

#' @export
#' @rdname indiv_inbreeding
indiv_inbreeding.tbl_df <- function(
    .x,
    method = c("WG17"),
    allele_sharing_mat = NULL,
    ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
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
    ...) {
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
indiv_inbreeding.grouped_df <- function(
    .x,
    method = c("WG17"),
    allele_sharing_mat = NULL,
    ...) {
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
