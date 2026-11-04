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
#'   "list", or "matrix". Default is "tidy".
#'   Data can be grouped by more than one variable: groups are then labelled
#'   by the values of the grouping variables pasted with "." (e.g. `pop1.F`),
#'   and the "tidy" output has one column per grouping variable rather than a
#'   single `group` column. The `group` column, or the columns of the
#'   grouping variables, keep the type of the grouping variables (e.g. a
#'   factor keeps its levels), so that groups stay in their order when the
#'   output is sorted or plotted.
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

  geno_fbm <- .gt_get_fbm(.x)
  # rows (individuals) that we want to use
  rows_to_keep <- .gt_fbm_rows(.x)

  # Keep the common diploid case on the ploidy-free kernel. Pseudohaploid
  # samples retain 0/2 dosage storage despite contributing one allele, so
  # they must use the ploidy-aware kernel. The cpp code reads the FBM
  # directly, in blocks of block_size loci, and fills a single matrix
  # (rather than binding the matrices from each block); openMP uses
  # multiple cores in the cpp code.
  if (is_diploid_only(.x)) {
    pi_mat <- gt_grouped_pi_diploid(
      BM = geno_fbm,
      rowInd = rows_to_keep,
      colInd = .gt_fbm_cols(.x),
      groupIds = dplyr::group_indices(.x) - 1L,
      ngroups = dplyr::n_groups(.x),
      ncores = n_cores,
      block_size = block_size
    )$pi
  } else {
    pi_mat <- gt_grouped_pi_polyploid(
      BM = geno_fbm,
      rowInd = rows_to_keep,
      colInd = .gt_fbm_cols(.x),
      groupIds = dplyr::group_indices(.x) - 1L,
      ngroups = dplyr::n_groups(.x),
      ploidy = indiv_ploidy(.x),
      is_pseudohaploid = is_pseudohaploid(.x),
      ncores = n_cores,
      block_size = block_size
    )$pi
  }

  groups <- .group_info(.x)
  pi_mat <- format_grouped_output(
    out_mat = pi_mat,
    group_ids = groups$labels,
    loci_names = loci_names(.x),
    type = type,
    group_keys = groups$keys
  )
  return(pi_mat)
}
