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
  if (operator == "custom") {
    if (is.null(f)) {
      stop("If operator is 'custom', f must be provided.")
    }
  }
  # check the loci table and the arguments that define the windows
  .check_window_args(
    loci_table,
    n_values = length(.x),
    window_size = window_size,
    step_size = step_size,
    size_unit = size_unit,
    min_loci = min_loci,
    complete = complete
  )
  # define the windows on each chromosome, and find the loci in each window
  # (this is done once, rather than scanning the whole loci table for each
  # chromosome)
  win <- .window_definitions(
    loci_table,
    window_size = window_size,
    step_size = step_size,
    size_unit = size_unit
  )
  # sum the values in each window (in the order of the loci, as runner does)
  win_sums <- window_sums_cpp(
    matrix(as.numeric(.x), ncol = 1),
    loci = win$loci,
    lo = win$lo,
    hi = win$hi,
    long_double = FALSE
  )
  win_sum <- win_sums$sum[, 1]
  win_n <- win_sums$n[, 1]

  if (operator != "custom") {
    # the mean (or sum) of the non-missing values; windows with no
    # values are NA (as in runner::mean_run() and runner::sum_run())
    if (operator == "mean") {
      stat <- win_sum / win_n
    } else {
      stat <- win_sum
    }
    stat[win_n == 0] <- NA
  } else {
    # apply the custom function with runner, one chromosome at a time
    stat <- vector("list", length(win$chrom_loci))
    for (i_chrom in seq_along(win$chrom_loci)) {
      chrom_win <- win$windows_chrom == i_chrom
      # note that runner measures the window from the end (not the front)
      stat[[i_chrom]] <- runner::runner(
        x = .x[win$chrom_loci[[i_chrom]]],
        k = window_size,
        at = win$windows$end[chrom_win],
        idx = if (size_unit == "bp") {
          loci_table$position[win$chrom_loci[[i_chrom]]]
        } else {
          integer(0)
        },
        na_pad = complete,
        f = f,
        ...
      )
    }
    stat <- unlist(stat)
  }
  # the number of loci with a non-missing value in each window (NA for
  # windows without loci)
  n_loci <- win_n
  n_loci[win$hi == win$lo] <- NA
  # if we only want complete windows, set incomplete ones to NA (as done by
  # na_pad in runner)
  if (complete) {
    stat[!win$complete] <- NA
    n_loci[!win$complete] <- NA
  }

  results <- win$windows
  results$stat <- stat
  results$n_loci <- n_loci
  # set values to NA if .n less then min_loci
  results$stat[results$n_loci < min_loci] <- NA

  return(results)
}

# Check the arguments that define windows (shared by the windows functions);
# n_values is the number of values (one per locus) to be summarised
.check_window_args <- function(loci_table, n_values, window_size, step_size,
                               size_unit, min_loci, complete) {
  # check that loci_table has the columns needed for this size_unit
  .check_loci_table_units(loci_table, size_unit)

  # check that complete is a boolean
  if (!is.logical(complete)) {
    stop("complete must be a boolean (logical).")
  }

  # check the loci_table has same number of rows as x
  if (nrow(loci_table) != n_values) {
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
  invisible(TRUE)
}

#' Define windows along each chromosome, and find the loci in each window
#'
#' Internal function used by the windows functions. On each chromosome (in the
#' order in which they appear in `loci_table`), windows are defined by their
#' end (`at`): the first window ends at the first multiple of `window_size`
#' that is not smaller than the index of the first locus, and windows then
#' move by `step_size` until the last locus is covered. As in `runner`, a
#' window ending at `at` includes the loci whose index (the position for
#' `size_unit = "bp"`, or the rank of the locus on its chromosome for
#' `size_unit = "snp"`) is in `(at - window_size, at]`.
#'
#' @param loci_table a loci table, with columns `chromosome` and (for
#'   `size_unit = "bp"`) `position`, sorted by position within each
#'   chromosome
#' @param window_size the size of the windows
#' @param step_size the step between windows
#' @param size_unit either "snp" or "bp"
#' @returns a list with elements:
#' - `windows`: a tibble with the `chromosome`, `start` and `end` of each
#'   window
#' - `windows_chrom`: the chromosome (as an index into `chrom_loci`) of each
#'   window
#' - `chrom_loci`: a list with the (1-based) indices of the loci of each
#'   chromosome, in the order of the loci table
#' - `loci`: the elements of `chrom_loci` concatenated
#' - `lo`, `hi`: window `w` includes the loci `loci[(lo[w] + 1):hi[w]]` (i.e.
#'   `lo` and `hi` are the 0-based, half-open range of the window in `loci`)
#' - `complete`: whether each window lies within the range of the loci of its
#'   chromosome (used to implement `complete = TRUE`)
#' - `locus_first`, `locus_last`: for each locus (in the order of the loci
#'   table), the first and last window that include it (0-based, with
#'   `locus_last < locus_first` if the locus is not in any window)
#' @keywords internal
#' @noRd
.window_definitions <- function(loci_table, window_size, step_size,
                                size_unit = c("snp", "bp")) {
  size_unit <- match.arg(size_unit)
  chrom <- loci_table$chromosome
  chromosomes <- unique(chrom)
  # the chromosome of each locus, as an index into chromosomes (for factors,
  # we match the integer codes, as matching factors would convert them to
  # character vectors)
  if (is.factor(chrom)) {
    chrom_id <- match(as.integer(chrom), as.integer(chromosomes))
  } else {
    chrom_id <- match(chrom, chromosomes)
  }
  # the indices of the loci of each chromosome, in the order of the table
  # (we split by a factor built directly from chrom_id, to avoid converting
  # it to character)
  chrom_factor <- structure(
    chrom_id,
    levels = as.character(seq_along(chromosomes)),
    class = "factor"
  )
  chrom_loci <- unname(split(seq_along(chrom), chrom_factor))
  # chromosome names as used in the results (factors are converted to
  # character, as when looping over them)
  if (is.factor(chromosomes)) {
    chromosomes <- as.character(chromosomes)
  }
  n_chrom <- length(chrom_loci)
  windows_at <- lo <- hi <- complete <- vector("list", n_chrom)
  locus_first <- locus_last <- integer(length(chrom))
  # number of loci and windows on the previous chromosomes
  loci_offset <- 0L
  win_offset <- 0L
  for (i_chrom in seq_len(n_chrom)) {
    if (size_unit == "bp") {
      position_sub <- loci_table$position[chrom_loci[[i_chrom]]]
    } else {
      position_sub <- seq_along(chrom_loci[[i_chrom]])
    }
    # window ends
    window_range <- ceiling(range(position_sub) / window_size)
    window_at <- seq(
      from = window_range[1] * window_size,
      to = window_range[2] * window_size,
      by = step_size
    )
    # if the last loci are beyond the last window (as the range is not a
    # multiple of step_size), add one more step to include them (unless they
    # fall in a gap between windows, when step_size is larger than
    # window_size)
    last_at <- window_at[length(window_at)]
    max_pos <- max(position_sub)
    if (last_at < max_pos && last_at + step_size - window_size < max_pos) {
      window_at <- c(window_at, last_at + step_size)
    }
    windows_at[[i_chrom]] <- window_at
    # runner works with integer indices, window ends and sizes
    idx <- as.integer(position_sub)
    at <- as.integer(window_at)
    k <- as.integer(window_size)
    if (is.unsorted(idx)) {
      stop("loci must be sorted by position within each chromosome")
    }
    # the number of loci on this chromosome before each window, and up to its
    # end: the loci in a window are those with idx > at - k and idx <= at
    lo_chrom <- findInterval(at - k, idx)
    hi_chrom <- findInterval(at, idx)
    lo[[i_chrom]] <- loci_offset + lo_chrom
    hi[[i_chrom]] <- loci_offset + hi_chrom
    # windows that extend beyond the first or last locus are incomplete
    complete[[i_chrom]] <- (at - k + 1 >= idx[1]) & (at <= idx[length(idx)])
    # for each locus (with rank p on the chromosome), the windows that
    # include it are those with lo_chrom < p <= hi_chrom; as lo_chrom and
    # hi_chrom are non-decreasing, these are a contiguous range of windows
    p <- seq_along(idx)
    locus_first[chrom_loci[[i_chrom]]] <-
      win_offset + findInterval(p - 1L, hi_chrom)
    locus_last[chrom_loci[[i_chrom]]] <-
      win_offset + findInterval(p - 1L, lo_chrom) - 1L
    loci_offset <- loci_offset + length(idx)
    win_offset <- win_offset + length(at)
  }
  n_windows <- lengths(windows_at)
  window_at <- unlist(windows_at)
  list(
    windows = tibble(
      chromosome = rep(chromosomes, n_windows),
      start = window_at - window_size + 1,
      end = window_at
    ),
    windows_chrom = rep(seq_len(n_chrom), n_windows),
    chrom_loci = chrom_loci,
    loci = unlist(chrom_loci),
    lo = unlist(lo),
    hi = unlist(hi),
    complete = unlist(complete),
    locus_first = locus_first,
    locus_last = locus_last
  )
}

# Check that a loci_table has the columns needed to define windows or blocks
# in a given size_unit ("snp" only needs 'chromosome', "bp" also needs
# 'position', and "cM" also needs a non-zero 'genetic_dist')
.check_loci_table_units <- function(loci_table, size_unit) {
  # if size_unit is snp, check that loci_table has column chromosome
  # else we also need position (or genetic_dist)
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
  } else if (size_unit == "cM") {
    if (!all(c("chromosome", "genetic_dist") %in% colnames(loci_table))) {
      stop(
        "loci_table must contain columns 'chromosome' and 'genetic_dist' ",
        "when size_unit is 'cM'."
      )
    }
    if (all(loci_table$genetic_dist == 0, na.rm = TRUE)) {
      stop(
        "all values of 'genetic_dist' are 0: genetic distances are needed ",
        "when size_unit is 'cM'."
      )
    }
    if (anyNA(loci_table$genetic_dist)) {
      stop("'genetic_dist' can not contain NAs when size_unit is 'cM'.")
    }
  }
  invisible(TRUE)
}
