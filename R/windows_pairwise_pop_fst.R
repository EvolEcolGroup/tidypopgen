#' Compute pairwise Fst for a sliding window
#'
#' @description This function computes pairwise Fst for a sliding window across
#'   each chromosome.
#'
#' @param .x a grouped `gen_tibble` object.
#'   If grouped by more than one variable, groups are labelled by the values
#'   of the grouping variables pasted with "." (e.g. `pop1.F`).
#' @param type type of object to return. One of "matrix" or "tidy". Default is
#'   "matrix". "matrix" returns a dataframe where each row is a window, followed
#'   by columns of Fst values for each pairwise population a and b comparison.
#'   "tidy" returns a tidy tibble of the same data in 'long' format, where each
#'   row is one window for one pairwise population a and b comparison.
#' @param method the method to use for calculating Fst, one of 'Hudson',
#'   'Nei87', and 'WC84'. See [pairwise_pop_fst()] for details.
#' @param window_size The size of the window to use for the estimates.
#' @param step_size The step size to use for the windows.
#' @param size_unit Either "snp" or "bp". If "snp", the window size and step
#'   size are in number of SNPs. If "bp", the window size and step size are in
#'   base pairs.
#' @param min_loci The minimum number of loci required to calculate a window
#'   statistic. If the number of loci in a window is less than this, the window
#'   statistic will be NA.
#' @param complete Should the function be evaluated on complete windows only? If
#'   FALSE, the default, then partial computations will be allowed at the end of
#'   the chromosome.
#' @returns either a data frame with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `start`: the starting locus of the window
#' - `end`: the ending locus of the window
#' - `fst_a.b`: the pairwise Fst value for the population a and b (there will be
#'   multiple such columns if there are more than two populations) or a tidy
#'   tibble with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `start`: the starting locus of the window
#' - `end`: the ending locus of the window
#' - `stat_name`: the name of population a and b used in the pairwise Fst
#'   calculation (e.g. "fst_pop1.pop2")
#' - `value`: the pairwise Fst value for the population a and b
#' @param ... other arguments passed to specific methods, currently
#'   unused.
#' @rdname windows_pairwise_pop_fst
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>%
#'   group_by(population) %>%
#'   windows_pairwise_pop_fst(
#'     window_size = 3, step_size = 2,
#'     size_unit = "snp", min_loci = 2
#'   )
#'
windows_pairwise_pop_fst <- function(.x, ...) {
  UseMethod("windows_pairwise_pop_fst", .x)
}

#' @export
#' @rdname windows_pairwise_pop_fst
windows_pairwise_pop_fst.tbl_df <- function(.x,
                                            type = c("matrix", "tidy"),
                                            method = c("Hudson", "Nei87", "WC84"),
                                            window_size,
                                            step_size,
                                            size_unit = c("snp", "bp"),
                                            min_loci = 1,
                                            complete = FALSE,
                                            ...) {
  rlang::check_dots_empty()
  # Check if the input is a gen_tibble
  stopifnot_gen_tibble(.x)
  method <- match.arg(method)
  type <- match.arg(type)
  size_unit <- match.arg(size_unit)

  # the same checks as in pairwise_pop_fst()
  if (!inherits(.x, "grouped_df")) {
    stop(".x should be a grouped df")
  }
  # only operate on diploids or, if pseudohaploids, Hudson
  .check_fst_ploidy(.x, method)
  # get the populations, and all pairs of populations
  groups <- .group_info(.x)
  pairwise_combn <- utils::combn(groups$n, 2)

  # define the windows, and find the windows that include each locus
  loci_table <- show_loci(.x)
  .check_window_args(
    loci_table,
    n_values = nrow(loci_table),
    window_size = window_size,
    step_size = step_size,
    size_unit = size_unit,
    min_loci = min_loci,
    complete = complete
  )
  win <- .window_definitions(
    loci_table,
    window_size = window_size,
    step_size = step_size,
    size_unit = size_unit
  )

  # sum the numerators and denominators of Fst of each pair of populations
  # in each window, in a single pass over the genotypes (without storing
  # them by locus)
  win_sums <- windows_pairwise_pop_fst_cpp(
    .gt_get_fbm(.x),
    rowInd = .gt_fbm_rows(.x),
    colInd = .gt_fbm_cols(.x),
    groupIds = groups$ids,
    ngroups = groups$n,
    ploidy = indiv_ploidy(.x),
    pairwise_combn = pairwise_combn,
    method = method,
    locus_first = win$locus_first,
    locus_last = win$locus_last,
    n_windows = nrow(win$windows),
    both_valid = FALSE,
    block_size = 1024,
    ncores = bigstatsr::nb_cores()
  )

  # mean of the numerators (or denominators) in each window, as computed by
  # windows_stats_generic() with operator = "mean": NA if the window has
  # fewer than min_loci values (or no values), or if it is incomplete and
  # only complete windows are requested
  window_means <- function(sums, n) {
    means <- sums / n
    means[n < min_loci | n == 0] <- NA
    if (complete) {
      means[!win$complete, ] <- NA
    }
    means
  }
  window_fst <- window_means(win_sums$sum_num, win_sums$n_num) /
    window_means(win_sums$sum_den, win_sums$n_den)

  # name the columns after the pairs of populations
  colnames(window_fst) <- paste0(
    "fst_",
    .pair_labels(groups$labels, pairwise_combn)
  )
  res <- dplyr::bind_cols(win$windows, tibble::as_tibble(window_fst))
  if (type == "matrix") {
    return(res)
  } else if (type == "tidy") {
    cols <- names(res)[!names(res) %in% c("chromosome", "start", "end")]
    res <-
      res %>% tidyr::pivot_longer(cols = all_of(cols), names_to = "stat_name")
    return(res)
  }
}
