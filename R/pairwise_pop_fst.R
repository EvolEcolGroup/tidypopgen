#' Compute pairwise population Fst
#'
#' This function computes pairwise Fst. The following methods are implemented:
#' - 'Hudson': Hudson's formulation, as derived in Bhatia et al (2013)
#' for diploids. This is the only method that can also be used with
#' pseudohaploid data.
#' - 'Nei87' : Fst according to Nei (1987) - includes the correction for
#' heterozygosity when computing Ht (it uses the same formulation as in
#' `hierfstat::pairwise.neifst()`),
#' - 'WC84' : Weir and Cockerham (1984), correcting for missing data (it uses
#' the same formulation as in `hierfstat::pairwise.WCfst()`).
#'
#' For all formulae, the genome wide estimate is obtained by taking the ratio of
#' the mean numerators and denominators over all relevant SNPs.
#'
#' @details Loci with no genotypes in one of the two populations of a pair
#'   (e.g. loci that failed in all the individuals of a population) carry no
#'   information on the differentiation between them, so they are dropped for
#'   that pair with all methods: their Fst by locus is `NaN`, and they do not
#'   contribute to the mean numerator and denominator of the genome wide
#'   estimate.
#'
#'   **Difference from `hierfstat` for 'WC84'**: [hierfstat::pairwise.WCfst()]
#'   (and [hierfstat::wc()]) keep these loci. Their among-population variance
#'   component (\eqn{\sigma^2_a}, the numerator) is set to 0, while their
#'   within-population components (\eqn{\sigma^2_b} and \eqn{\sigma^2_w},
#'   estimated from the one population with data) still add to the
#'   denominator. So `hierfstat` gives these loci an Fst of 0 rather than
#'   `NaN`, and its genome wide estimate is pulled towards 0 (i.e. it is
#'   usually lower than the one computed here).
#'   The two estimates are identical when every locus has data in both
#'   populations, or if such loci are removed before using `hierfstat`. This
#'   applies to every function that computes Fst with 'WC84' (e.g.
#'   [windows_pairwise_pop_fst()], [jackknife_pairwise_pop_fst()] and
#'   [nwise_pop_pbs()]). 'Nei87' gives the same results as
#'   [hierfstat::pairwise.neifst()], which also drops these loci.
#'
#' @references Bhatia G, Patterson N, Sankararaman S, Price AL. (2013)
#'   Estimating and Interpreting FST: The Impact of Rare Variants. Genome
#'   Research, 23(9):1514–1521.
#'
#'   Nei, M. (1987) Molecular Evolutionary Genetics. Columbia University Press
#'
#'   Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the
#'   analysis of population structure. Evolution, 38(6): 1358–1370.
#'
#' @param .x a grouped [`gen_tibble`] (as obtained by using [dplyr::group_by()])
#' @param type type of object to return One of "tidy" or "pairwise" for a
#'   pairwise matrix of populations. Default is "tidy". If data are grouped by
#'   more than one variable, the "tidy" output has one column per grouping
#'   variable for each population of a pair (e.g. `population_1`, `sex_1`,
#'   `population_2`, `sex_2`), and populations are otherwise labelled by the
#'   values of the grouping variables pasted with "." (e.g. `pop1.F`).
#' @param by_locus_type type of object to return. One of "tidy", "matrix" or
#'   "list". Default is "tidy".
#' @param by_locus boolean, determining whether Fst should be returned by
#'   locus(TRUE), or as a single genome wide value obtained by taking the ratio
#'   of the mean numerator and denominator (FALSE, the default).
#' @param method one of 'Hudson', 'Nei87', and 'WC84'
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param return_num_dem returns a list of numerators and denominators for each
#'   locus. This is useful for creating windowed estimates of Fst (as we need to
#'   compute the mean numerator and denominator within each window). Default is
#'   FALSE.
#' @returns if `type=tidy`, a tibble of genome-wide pairwise Fst values with
#'   each pairwise combination as a row if "by_locus=FALSE", else a list
#'   including the tibble of genome-wide values as well as a matrix with
#'   pairwise Fst by locus with loci as rows and and pairwise combinations as
#'   columns. If `type=pairwise`, a matrix of genome-wide pairwise Fst values is
#'   returned.
#' @param ... other arguments passed to specific methods, currently
#'   unused.
#' @rdname pairwise_pop_fst
#' @export
#' @seealso [hierfstat::pairwise.neifst()]
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # For a basic global pairwise Fst calculation:
#' example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87")
#'
#' # With a pairwise matrix:
#' example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87", type = "pairwise")
#'
#' # To calculate Fst by locus:
#' example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Hudson", by_locus = TRUE)
#'
pairwise_pop_fst <- function(.x, ...) {
  UseMethod("pairwise_pop_fst", .x)
}

#' @export
#' @rdname pairwise_pop_fst
pairwise_pop_fst.gen_tbl <- function(
  .x,
  type = c("tidy", "pairwise"),
  by_locus = FALSE,
  by_locus_type = c("tidy", "matrix", "list"),
  method = c("Hudson", "Nei87", "WC84"),
  return_num_dem = FALSE,
  n_cores = bigstatsr::nb_cores(),
  ...
) {
  rlang::check_dots_empty()
  if (n_cores > 1) {
    # Remove checking for two levels of parallelism
    .old_opt <- getOption("bigstatsr.check.parallel.blas", TRUE)
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = .old_opt), add = TRUE)
  }

  type <- match.arg(type)
  if (by_locus) {
    by_locus_type <- match.arg(by_locus_type)
  }

  if (!inherits(.x, "grouped_df")) {
    stop(".x should be a grouped df")
  }
  method <- match.arg(method)

  if (!is.logical(return_num_dem)) {
    stop("return_num_dem must be a logical value (TRUE or FALSE)")
  }

  # if we want to return the numerator and denominator, we need to
  # compute the pairwise Fst for each locus
  if (return_num_dem && !by_locus) {
    message("`by_locus` set to TRUE because `return_num_dem = TRUE`.")
    by_locus <- TRUE
  }

  # only operate on diploids or, if pseudohaploids, Hudson
  .check_fst_ploidy(.x, method)

  # get the populations
  groups <- .group_info(.x)
  # create all combinations
  pairwise_combn <- utils::combn(groups$n, 2)

  # compute the numerator and denominator of Fst for every pair of
  # populations in a single pass over the genotypes: loci are processed in
  # blocks, summarising the populations within each block and then looping
  # over the pairs, so that we never store the full loci x populations
  # summaries (only the matrices by locus are stored, if requested)
  fst_list <- pairwise_pop_fst_cpp(
    .gt_get_fbm(.x),
    rowInd = .gt_fbm_rows(.x),
    colInd = .gt_fbm_cols(.x),
    groupIds = groups$ids,
    ngroups = groups$n,
    ploidy = indiv_ploidy(.x),
    pairwise_combn = pairwise_combn,
    method = method,
    by_locus = by_locus,
    return_num_dem = return_num_dem,
    block_size = 1024,
    ncores = n_cores
  )

  # the groups of each pair, with one column per grouping variable for each
  # of the two groups (e.g. population_1 and population_2)
  group_combinations <- cbind(
    groups$keys[pairwise_combn[1, ], ],
    groups$keys[pairwise_combn[2, ], ]
  )
  names(group_combinations) <- c(
    paste0(groups$vars, "_1"),
    paste0(groups$vars, "_2")
  )
  # the names of the pairs (e.g. pop1.pop2), used to name columns
  pair_names <- .pair_labels(groups$labels, pairwise_combn)
  # if we want the numerator and denominator,
  # we need to format them and return them
  if (return_num_dem) {
    # TODO turn this into matrices with numerical row names
    fst_list$Fst_by_locus_num <- matrix_int_names(
      data = fst_list$Fst_by_locus_num,
      row_names = loci_names(.x),
      col_names = pair_names
    )
    fst_list$Fst_by_locus_den <- matrix_int_names(
      data = fst_list$Fst_by_locus_den,
      row_names = loci_names(.x),
      col_names = pair_names
    )

    return(fst_list)
  }

  # otherwise we start formatting the other objects
  fst_tot <- tibble::tibble(group_combinations, value = fst_list$fst_tot)
  # now add a class
  class(fst_tot) <- c("pairwise_tbl", class(fst_tot))

  if (type == "pairwise") { # if we return a matrix
    fst_tot <- .pairwise_matrix(fst_list$fst_tot, groups$labels, pairwise_combn)
  }

  if (by_locus && by_locus_type == "matrix") {
    # turn matrix into matrix_int
    fst_list$fst_locus <- matrix_int_names(
      data = fst_list$fst_locus,
      row_names = loci_names(.x),
      col_names = paste0("fst_", pair_names)
    )

    return(list(Fst_by_locus = fst_list$fst_locus, Fst = fst_tot))
  } else if (by_locus && by_locus_type == "tidy") {
    fst_mat_tbl <- as.data.frame(fst_list$fst_locus)
    colnames(fst_mat_tbl) <- paste0("fst_", pair_names)
    fst_mat_tbl$loci <- loci_names(.x)
    cols <- names(fst_mat_tbl)[names(fst_mat_tbl) != "loci"]
    long_fst_loc <- fst_mat_tbl %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(cols),
        names_to = "stat_name"
      )
    return(list(Fst_by_locus = long_fst_loc, Fst = fst_tot))
  } else if (by_locus && by_locus_type == "list") {
    fst_mat_tbl <- as.data.frame(fst_list$fst_locus)
    colnames(fst_mat_tbl) <- paste0("fst_", pair_names)
    fst_list$fst_locus <- as.list(fst_mat_tbl)
    return(list(Fst_by_locus = fst_list$fst_locus, Fst = fst_tot))
  } else {
    return(fst_tot)
  }
}

# Turn the values for each pair of groups into a symmetric matrix (with NA
# on the diagonal), with the group labels as row and column names
.pairwise_matrix <- function(values, labels, pairwise_combn) {
  n <- length(labels)
  mat <- matrix(NA_real_, n, n, dimnames = list(labels, labels))
  mat[t(pairwise_combn)] <- values
  mat[t(pairwise_combn[2:1, , drop = FALSE])] <- values
  class(mat) <- c("pairwise_matrix", class(mat))
  mat
}

# Check that Fst can be computed on the ploidy of .x with the given method:
# Fst is only implemented for diploids or, if the data include
# pseudohaploids, for Hudson (used by pairwise_pop_fst(),
# windows_pairwise_pop_fst() and the jackknife functions)
.check_fst_ploidy <- function(.x, method) {
  if (!is_diploid_only(.x) && !is_pseudohaploid(.x)) {
    stop("fst can be currently computed on diploid data")
  }
  if (is_pseudohaploid(.x) && method != "Hudson") {
    stop("only `method = Hudson` is valid when the data include pseudohaploids")
  }
  invisible(TRUE)
}

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
pairwise_pop_fst.data.frame <- function(.x, ...) {
  pairwise_pop_fst(restore_genotypes_tbl(.x, "pairwise_pop_fst"), ...)
}
