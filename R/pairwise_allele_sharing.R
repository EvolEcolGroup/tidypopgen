#' Compute the Pairwise Allele Sharing Matrix for a `gen_tibble` object
#'
#' This function computes the Allele Sharing matrix. Estimates Allele Sharing
#' (equivalent to the quantity estimated by `hierfstat::matching()`) between
#' pairs of individuals (for each locus, gives 1 if the two individuals are
#' homozygous for the same allele, 0 if they are homozygous for a different
#' allele, and 1/2 if at least one individual is heterozygous. Matching is the
#' average of these 0, 1/2 and 1s)
#'
#' Imputed genotypes are used if `x` is set to use them (see
#' [gt_set_imputed()]), and treated as missing otherwise.
#'
#' `pairwise_allele_sharing()` is an S3 generic. The method for a
#' [`gen_tibble`] computes the matrix with the method for its `genotypes`
#' column (of class `vctrs_bigSNP`), and labels it with the `id` of the
#' individuals; grouping is ignored. Other packages can add methods for their
#' own classes of genotypes (e.g. genotype likelihoods).
#'
#' @param x a `gen_tibble` object, or a vector of class `vctrs_bigSNP` (the
#'   `genotypes` column of a [`gen_tibble`], whose individuals are then
#'   labelled with the names of the vector).
#' @param as_matrix boolean, determining whether the results should be a square
#'   symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)
#' @param block_size maximum number of loci read at once. Genotypes are
#'   compared in a compact bit-packed form, so this has little effect on
#'   memory use.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns a matrix of allele sharing between all pairs of individuals (the
#'   diagonal is the allele sharing of each individual with itself, 1 minus
#'   half its proportion of heterozygous loci), or a tidy tibble
#' @rdname pairwise_allele_sharing
#' @export
#' @seealso [hierfstat::matching()]
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Compute allele sharing between individuals
#' example_gt %>% pairwise_allele_sharing(as_matrix = FALSE)
#'
#' # Alternatively, return as a tibble
#' example_gt %>% pairwise_allele_sharing(as_matrix = TRUE)
pairwise_allele_sharing <- function(x, ...) {
  UseMethod("pairwise_allele_sharing", x)
}

#' @export
#' @rdname pairwise_allele_sharing
pairwise_allele_sharing.tbl_df <- function(
  x,
  as_matrix = FALSE,
  block_size = bigstatsr::block_size(nrow(x)),
  n_cores = bigstatsr::nb_cores(),
  ...
) {
  stopifnot_has_genotypes(x)
  # hand the genotypes column over to its method, which returns a
  # pairwise_matrix, and label it with the ids of the tibble (the matrix is
  # modified in place, without copies)
  ashare_matrix <- pairwise_allele_sharing(
    x$genotypes,
    as_matrix = TRUE,
    block_size = block_size,
    n_cores = n_cores,
    ...
  )
  dimnames(ashare_matrix) <- list(x$id, x$id)
  if (as_matrix) {
    return(ashare_matrix)
  } else {
    class(ashare_matrix) <- NULL
    return(tidy_dist_matrix(ashare_matrix))
  }
}

#' @export
#' @rdname pairwise_allele_sharing
pairwise_allele_sharing.vctrs_bigSNP <- function(
  x,
  as_matrix = FALSE,
  block_size = bigstatsr::block_size(length(x)),
  n_cores = bigstatsr::nb_cores(),
  ...
) {
  rlang::check_dots_empty()
  stopifnot_diploid(x)
  ashare_matrix <- snp_allele_sharing(
    attr(x, "fbm"),
    ind.row = vctrs::vec_data(x),
    ind.col = show_loci(x)$big_index,
    block.size = block_size,
    n_cores = n_cores
  )
  # the individuals are labelled with the names of the vector (the ids of the
  # gen_tibble)
  dimnames(ashare_matrix) <- list(names(x), names(x))
  if (as_matrix) {
    # add a class to indicate this is a pairwise matrix
    class(ashare_matrix) <- c("pairwise_matrix", class(ashare_matrix))
    return(ashare_matrix)
  } else {
    return(tidy_dist_matrix(ashare_matrix))
  }
}
