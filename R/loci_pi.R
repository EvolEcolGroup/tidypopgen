#' Estimate nucleotide diversity (pi) at each locus
#'
#' Estimate nucleotide diversity (\eqn{\pi}) at each locus, accounting for
#' missing values.
#'
#' Nucleotide diversity is estimated as
#'
#' \deqn{
#' \pi = \frac{c_0 c_1}{n(n - 1)/2},
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{c_0} is the number of observed copies of allele 0,
#'   \item \eqn{c_1} is the number of observed copies of allele 1,
#'   \item \eqn{n = c_0 + c_1} is the total number of non-missing allele
#'     copies at the locus.
#' }
#'
#' Missing values
#' are accommodated by computing \eqn{n}, \eqn{c_0}, and \eqn{c_1} using only
#' non-missing genotype data. For polyploid data, each individual contributes
#' a number of allele copies equal to its ploidy.
#'
#' This corresponds to the unbiased gene diversity estimator of
#' Nei and Roychoudhury (1974) and is applicable to diploid, polyploid
#' mixed-ploidy, and pseudohaploid `gen_tibble` objects alike.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#' directly to the function). This defaults to "genotypes" and can only take
#' that value. There is no need for the user to set it, but it is included to
#' resolve certain tidyselect operations.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param block_size maximum number of loci read at once.
#' @param type type of object to return, if using grouped method. One of "tidy",
#' "list", or "matrix". Default is "tidy".
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_pi
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # For pi
#' example_gt %>% loci_pi()
#'
#' # For pi per locus per population
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_pi()
#' # alternatively, return a list of populations with their pi
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_pi(type = "list")
#' # or a matrix with populations in columns and loci in rows
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_pi(type = "matrix")
#' # or within reframe (not recommended, as it much less efficient
#' # than using it directly as shown above)
#' example_gt %>%
#'   group_by(population) %>%
#'   reframe(pi = loci_pi(genotypes))
loci_pi <- function(.x, .col = "genotypes", n_cores, block_size, type, ...) {
  UseMethod("loci_pi", .x)
}

#' @export
#' @rdname loci_pi
loci_pi.tbl_df <- function(
  .x,
  .col = "genotypes",
  # multicore is used by openMP within the
  # freq cpp function
  n_cores = bigstatsr::nb_cores(),
  # the bigapply that splits in blocks is not
  # multithreaded, as we use the multiple threads
  # for openMP
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
) {
  # TODO this is a hack to deal with the class being dropped when going
  # through group_map
  stopifnot_gen_tibble(.x)
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_pi only works with the genotypes column")
  }
  loci_pi(.x$genotypes, n_cores = n_cores, block_size = block_size)
}


#' @export
#' @rdname loci_pi
loci_pi.vctrs_bigSNP <- function(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
) {
  rlang::check_dots_empty()

  # get the FBM
  geno_fbm <- attr(.x, "fbm")
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep) > 1) {
    # Keep the common diploid case on the ploidy-free kernel. Pseudohaploid
    # samples retain 0/2 dosage storage despite contributing one allele, so
    # they must use the ploidy-aware kernel. The cpp code reads the FBM
    # directly, in blocks of block_size loci, using openMP for
    # parallelisation.
    if (is_diploid_only(.x)) {
      pi <- gt_pi_diploid(
        BM = geno_fbm,
        rowInd = rows_to_keep,
        colInd = attr(.x, "loci")$big_index,
        ncores = n_cores,
        block_size = block_size
      )
    } else {
      pi <- gt_pi_polyploid(
        BM = geno_fbm,
        rowInd = rows_to_keep,
        colInd = attr(.x, "loci")$big_index,
        ploidy = indiv_ploidy(.x),
        is_pseudohaploid = is_pseudohaploid(.x),
        ncores = n_cores,
        block_size = block_size
      )
    }
  } else {
    # if we have a single individual
    # pi does not really make sense for a single individual
    pi <- rep(NA_real_, length(attr(.x, "loci")$big_index))
  }
  pi
}

#' @export
#' @rdname loci_pi
loci_pi.grouped_df <- function(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  type = c("tidy", "list", "matrix"),
  ...
) {
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_pi only works with the genotypes column")
  }
  rlang::check_dots_empty()
  type <- match.arg(type)
  # if we only have one individual, return NA for all loci
  if (nrow(.x) == 1) {
    return(rep(NA_real_, nrow(show_loci(.x))))
  }

  # check that we only have one grouping variable
  if (length(.x %>% dplyr::group_vars()) > 1) {
    stop("loci_pi only works with one grouping variable")
  }

  # compute pi for each group (the number of valid allele copies is not
  # needed, so it is not returned)
  pi_mat <- .loci_pi_n(
    .x$genotypes,
    group_ids = dplyr::group_indices(.x) - 1L,
    n_groups = dplyr::n_groups(.x),
    n_cores = n_cores,
    block_size = block_size,
    return_n = FALSE
  )$pi

  pi_mat <- format_grouped_output(
    out_mat = pi_mat,
    group_ids = dplyr::group_keys(.x) %>% pull(1),
    loci_names = loci_names(.x),
    type = type
  )
  return(pi_mat)
}

#' Estimate pi and the number of valid allele copies at each locus by group
#'
#' Internal function that computes nucleotide diversity (pi) at each locus for
#' each group of individuals, together with the number of non-missing allele
#' copies (`n`) used to estimate pi at that locus. Pi is `NA` when a group has
#' fewer than 2 valid allele copies at a locus. Ungrouped data can be handled
#' by assigning all individuals to a single group (i.e. `group_ids` set to 0
#' for all individuals, and `n_groups = 1`).
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column
#'   of a [`gen_tibble`] object).
#' @param group_ids a vector of 0-based group indices, one per individual in
#'   `.x`.
#' @param n_groups the number of groups.
#' @param n_cores number of cores to be used.
#' @param block_size maximum number of loci read at once.
#' @param return_n boolean, should the number of valid allele copies be
#'   returned? If FALSE, `n` is NULL (the matrix is as large as `pi`, so we
#'   only build it when needed).
#' @returns a list with two matrices with one row per locus and one column per
#'   group: `pi`, the nucleotide diversity, and `n`, the number of valid
#'   allele copies.
#' @keywords internal
#' @noRd
.loci_pi_n <- function(.x, group_ids, n_groups, n_cores, block_size,
                       return_n = TRUE) {
  geno_fbm <- attr(.x, "fbm")
  # rows (individuals) and columns (loci) of the FBM that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  cols_to_keep <- attr(.x, "loci")$big_index

  # Keep the common diploid case on the ploidy-free kernel. Pseudohaploid
  # samples retain 0/2 dosage storage despite contributing one allele, so
  # they must use the ploidy-aware kernel. The cpp code reads the FBM
  # directly, in blocks of block_size loci, and fills a single matrix
  # (rather than binding the matrices from each block); openMP uses
  # multiple cores in the cpp code.
  if (is_diploid_only(.x)) {
    pi_n <- gt_grouped_pi_diploid(
      BM = geno_fbm,
      rowInd = rows_to_keep,
      colInd = cols_to_keep,
      groupIds = group_ids,
      ngroups = n_groups,
      ncores = n_cores,
      block_size = block_size,
      return_n = return_n
    )
  } else {
    pi_n <- gt_grouped_pi_polyploid(
      BM = geno_fbm,
      rowInd = rows_to_keep,
      colInd = cols_to_keep,
      groupIds = group_ids,
      ngroups = n_groups,
      ploidy = indiv_ploidy(.x),
      is_pseudohaploid = is_pseudohaploid(.x),
      ncores = n_cores,
      block_size = block_size,
      return_n = return_n
    )
  }
  list(pi = pi_n$pi, n = pi_n$n)
}
