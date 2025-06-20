#' Estimate window statistics from per locus estimates
#'
#' This function is mostly designed for developers: it is a general function to
#' estimate window statistics from per locus estimates. This function takes a
#' vector of per locus estimates, and aggregates them by sum or mean per window.
#' To compute specific quantities directly from a `gen_tibble`, use the
#' appropriate `window_*` functions, e.g [windows_pairwise_pop_fst()] to compute
#' pairwise Fst.
#'
#' @param .x A vector containing the per locus estimates.
#' @param loci_table a dataframe including at least a column 'chromosome', and
#'   additionally a column 'position' if `size_unit` is "bp".
#' @param operator The operator to use for the window statistics. Either "mean",
#'   "sum" or "custom" to use a custom function `.f`.
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
#' @param f a custom function to use for the window statistics. This function
#'   should take a vector of locus estimates and return a single value.
#' @param ... Additional arguments to be passed to the custom operator function.
#' @returns A tibble with columns: 'chromosome', 'start', 'end', 'stats', and
#'   'n_loci'. The 'stats' column contains the mean of the per locus estimates
#'   in the window, and 'n_loci' contains the number of loci in the window.
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' miss_by_locus <- loci_missingness(example_gt)
#'
#' # Calculate mean missingness across windows
#' windows_stats_generic(miss_by_locus,
#'   loci_table = show_loci(example_gt),
#'   operator = "mean", window_size = 1000,
#'   step_size = 1000, size_unit = "bp",
#'   min_loci = 1, complete = FALSE
#' )
#'
windows_stats_generic <- function(.x, loci_table,
                                  operator = c("mean", "sum", "custom"),
                                  window_size, step_size,
                                  size_unit = c("snp", "bp"), min_loci = 1,
                                  complete = FALSE,
                                  f = NULL, ...) {
  size_unit <- match.arg(size_unit)
  operator <- match.arg(operator)
  if (operator == "mean") {
    operator_func <- runner::mean_run
  } else if (operator == "sum") {
    operator_func <- runner::sum_run
  } else if (operator == "custom") {
    if (is.null(f)) {
      stop("If operator is 'custom', f must be provided.")
    }
    operator_func <- runner::runner
  }
  # if size_unit is snp, check that loci_table has column chromosome
  # else we also need position
  if (size_unit == "snp") {
    if (!all(c("chromosome") %in% colnames(loci_table))) {
      stop(
        "loci_table must contain column 'chromosome' when ",
        "size_unit is 'snp'."
      )
    }
  } else if (size_unit == "bp") {
    if (!all(c("chromosome", "position") %in% colnames(loci_table))) {
      stop(
        "loci_table must contain columns 'chromosome' and 'position' ",
        "when size_unit is 'bp'."
      )
    }
  }

  # check that complete is a boolean
  if (!is.logical(complete)) {
    stop("complete must be a boolean (logical).")
  }


  # check the loci_table has same number of rows as x
  if (nrow(loci_table) != length(.x)) {
    stop("loci_table must have the same number of rows as x.")
  }

  # check that window_size and step_size are positive
  if (window_size <= 0) {
    stop("window_size must be positive.")
  }
  if (step_size <= 0) {
    stop("step_size must be positive.")
  }
  # check that min_loci is positive
  if (min_loci <= 0) {
    stop("min_loci must be positive.")
  }
  # check that min_loci is less than window_size
  if (min_loci > window_size) {
    stop("min_loci must be less than window_size.")
  }
  # get unique chromosomes
  chromosomes <- unique(loci_table$chromosome)
  # create empty object to store results
  results <- NULL
  # loop over chromosomes
  for (i_chrom in chromosomes) {
    # get the stats for this chromosome
    x_sub <- .x[loci_table$chromosome == i_chrom]
    if (size_unit == "bp") {
      position_sub <- loci_table$position[loci_table$chromosome == i_chrom]
      idx_pos <- position_sub # use the positions as indices
    } else {
      position_sub <- seq_along(x_sub)
      idx_pos <- integer(0) # no need for indices, we just count the snps
    }
    # window ends (ranger uses the end as its `at` parameter)
    window_range <- ceiling(range(position_sub) / window_size)
    window_at <- seq(
      from = window_range[1] * window_size,
      to = window_range[2] * window_size,
      by = step_size
    )

    # note that runner measures the window from the end (not the front)
    if (operator != "custom") {
      stat <- operator_func(
        x = x_sub,
        k = window_size,
        at = window_at,
        idx = idx_pos,
        na_rm = TRUE,
        na_pad = complete
      )
    } else {
      stat <- operator_func(
        x = x_sub,
        k = window_size,
        at = window_at,
        idx = idx_pos,
        na_pad = complete,
        f = f,
        ...
      )
    }


    res <- tibble(
      chromosome = i_chrom,
      start = window_at - window_size + 1,
      end = window_at,
      stat = stat,
      n_loci = runner::sum_run(
        x = !is.na(x_sub),
        k = window_size,
        at = window_at,
        idx = idx_pos,
        na_rm = TRUE,
        na_pad = complete
      )
    )

    if (is.null(results)) { # for the first chromosome
      results <- res
    } else {
      results <- dplyr::bind_rows(results, res)
    }
  }
  # set values to NA if .n less then min_loci
  results$stat[results$n_loci < min_loci] <- NA

  return(results)
}
