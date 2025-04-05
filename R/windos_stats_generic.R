#' Estimate window statistics from per locus estimates
#'
#' This is a general function to estimate window statistics from per locus
#' estimates.
#'
#' @param x A vector containing the per locus estimates.
#' @param loci_table a dataframe with columns: 'chromosome', 'position', and
#'   'name'.
#' @param window_size The size of the window to use for the estimates.
#' @param window_step The step size to use for the windows.
#' @param size_unit Either "snp" or "bp". If "snp", the window size and step
#'   size are in number of SNPs. If "bp", the window size and step size are in
#'   base pairs.
#' @param min_loci The minimum number of loci required to calculate a window
#'   statistic. If the number of loci in a window is less than this, the window
#'   statistic will be NA.
#' @param complete Should the function be evaluated on complete windows only?
#' If FALSE, the default, then partial computations will be allowed at the end
#' of the chromosome.
#' @returns A dataframe with columns: 'chromosome', 'start', 'end', 'stats', and
#'   'n_snps'. The 'stats' column contains the mean of the per locus estimates
#'   in the window, and 'n_snps' contains the number of loci in the window.
#' @export

window_stats_generic <- function(x, loci_table, operator = c("mean", "sum"),
                                 window_size, window_step,
                                 size_unit = c("snp", "bp"), min_loci = 1,
                                 complete = FALSE) {
  size_unit <- match.arg(size_unit)
  operator <- match.arg(operator)
  operator_func <- slider_operator(operator, size_unit)
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

  # check the loci_table has same number of rows as x
  if (nrow(loci_table) != length(x)) {
    stop("loci_table must have the same number of rows as x.")
  }

  # check that window_size and window_step are positive
  if (window_size <= 0) {
    stop("window_size must be positive.")
  }
  if (window_step <= 0) {
    stop("window_step must be positive.")
  }
  # check that min_loci is positive
  if (min_loci <= 0) {
    stop("min_loci must be positive.")
  }
  # check that min_loci is less than window_size
  if (min_loci >= window_size) {
    stop("min_loci must be less than window_size.")
  }
  # loop over chromosomes
  # get unique chromosomes
  chromosomes <- unique(loci_table$chromosome)
  # create empty object to store results
  results <- NULL
  # loop over chromosomes
  for (i_chrom in chromosomes) {
    # get the stats for this chromosome
    x_sub <- x[loci_table$chromosome == i_chrom]
    if (size_unit == "snp") {
      res <- chrom_window_by_snp(
        x_sub = x_sub,
        i_chrom = i_chrom,
        operator_func = operator_func,
        window_size = window_size,
        window_step = window_step,
        min_loci = min_loci,
        complete = complete
      )
    } else if (size_unit == "bp") {
      position_sub <- loci_table$position[loci_table$chromosome == i_chrom]
      stop("windows with size unit bp have not been implemented yet")
      res <- chrom_window_by_pos(
        x_sub = x_sub,
        position_sub = position_sub,
        i_chrom = i_chrom,
        operator_func = operator_func,
        window_size = window_size,
        window_step = window_step,
        min_loci = min_loci,
        complete = complete
      )
    }

    if (is.null(results)) { # for the first chromosome
      results <- res
    } else {
      results <- rbind(results, res)
    }
  }
  return(results)
}


#' Internal function to choose the slider operator for windowing
#'
#' @param operator The operator to use for the window statistics.
#' @param size_unit The unit of the window size. Either "snp" or "bp".
#' @returns The function to use for the window statistics.
#' @keywords internal

slider_operator <- function(operator, size_unit) {
  if (operator == "mean") {
    if (size_unit == "snp") {
      operator_func <- slider::slide_mean
    } else if (size_unit == "bp") {
      operator_func <- slider::slide_index_mean
    }
  } else if (operator == "sum") {
    if (size_unit == "snp") {
      operator_func <- slider::slide_sum
    } else if (size_unit == "bp") {
      operator_func <- slider::slide_index_sum
    }
  }
  return(operator_func)
}



#' Internal function to compute the windows over a chromosome by snps
#' @param x_sub A vector containing the per locus estimates for a single
#' chromosome.
#' @param i_chrom The chromosome label
#'
#' @param operator_func The function to use for the window statistics. This
#' should be a function from the slider package, such as slide_mean or
#' slide_sum.
#' @param window_size The size of the window to use for the estimates.
#' @param window_step The step size to use for the windows.
#' @param min_loci The minimum number of loci required to calculate a window
#' statistic. If the number of loci in a window is less than this, the window
#' statistic will be NA.
#' @param complete Should the function be evaluated on complete windows only?
#' If FALSE, the default, then partial computations will be allowed at the end
#' of the chromosome.
#' @returns A dataframe with columns: 'chromosome', 'start', 'end', 'stats', and
#' 'n_snps'. The 'stats' column contains the mean of the per locus estimates
#' in the window, and 'n_snps' contains the number of loci in the window.
#' @keywords internal

chrom_window_by_snp <- function(x_sub, i_chrom, operator_func,
                                window_size, window_step,
                                min_loci, complete) {
  # get positions min and max for ranges
  this_ranges <- slider::slide(seq_along(x_sub),
    .f = range,
    .after = window_size - 1,
    .step = window_step,
    .complete = complete
  )
  # only some snps are starting point for a window
  real_windows <- unlist(lapply(this_ranges, function(x) !is.null(x)))
  # use slider to do window operations
  mean_values <- operator_func(x_sub,
    after = window_size - 1,
    step = window_step,
    complete = complete,
    na_rm = TRUE
  )
  n_valid <- slider::slide_sum(!is.na(x_sub),
    after = window_size - 1,
    step = window_step,
    complete = complete,
    na_rm = TRUE
  )
  mean_values[n_valid < min_loci] <- NA
  # subset to only the existing windows
  this_ranges <- this_ranges[real_windows]
  this_ranges <- do.call(rbind, this_ranges)
  res <- tibble::tibble(
    chromosome = i_chrom,
    start = this_ranges[, 1],
    end = this_ranges[, 2],
    stats = mean_values[real_windows],
    n_snps = n_valid[real_windows]
  )
  return(res)
}


#' Internal function to compute the windows over a chromosome by position
#' @param x_sub A vector containing the per locus estimates for a single
#' chromosome.
#' @param i_chrom The chromosome label
#' @param operator_func The function to use for the window statistics. This
#' should be a function from the slider package, such as slide_mean or
#' slide_sum.
#' @param window_size The size of the window to use for the estimates.
#' @param window_step The step size to use for the windows.
#' @param min_loci The minimum number of loci required to calculate a window
#' statistic. If the number of loci in a window is less than this, the window
#' statistic will be NA.
#' @param complete Should the function be evaluated on complete windows only?
#' If FALSE, the default, then partial computations will be allowed at the end
#' of the chromosome.
#' @returns A dataframe with columns: 'chromosome', 'start', 'end', 'stats', and
#' 'n_snps'. The 'stats' column contains the mean of the per locus estimates
#' in the window, and 'n_snps' contains the number of loci in the window.
#' @keywords internal

chrom_window_by_pos <- function(x_sub, position_sub, i_chrom,
                                operator_func,
                                window_size, window_step,
                                min_loci, complete) {
  # get positions min and max for ranges
  this_ranges <- slider::slide_index(position_sub,
    .i = position_sub,
    .f = range,
    .after = window_size - 1,
    .step = window_step,
    .complete = complete
  )
  # only some snps are starting point for a window
  real_windows <- unlist(lapply(this_ranges, function(x) !is.null(x)))
  # use slider to do window operations
  browser()
  mean_values <- operator_func(x_sub,
    i = position_sub,
    after = window_size - 1,
    step = window_step,
    complete = complete,
    na_rm = TRUE
  )
  n_valid <- slider::slide_index_sum(!is.na(x_sub),
    i = position_sub,
    after = window_size - 1,
    step = window_step,
    complete = complete,
    na_rm = TRUE
  )
  mean_values[n_valid < min_loci] <- NA
  # subset to only the existing windows
  this_ranges <- this_ranges[real_windows]
  this_ranges <- do.call(rbind, this_ranges)
  res <- tibble::tibble(
    chromosome = i_chrom,
    start = this_ranges[, 1],
    end = this_ranges[, 2],
    stats = mean_values[real_windows],
    n_snps = n_valid[real_windows]
  )
  return(res)
}
