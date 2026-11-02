#' Estimate a block jackknife standard error from per locus estimates
#'
#' This function is mostly designed for developers: it is a general function to
#' estimate the genome-wide value of a statistic, together with its standard
#' error, using a weighted block jackknife. It takes a vector of per locus
#' estimates, and aggregates them by sum or mean (or a custom function),
#' leaving out one block of loci at a time. To compute specific quantities
#' directly from a `gen_tibble`, use the appropriate `jackknife_*` functions,
#' e.g [jackknife_pairwise_pop_fst()] to compute pairwise Fst.
#'
#' @details The genome is divided into non-overlapping blocks, which never span
#'   more than one chromosome. Blocks can be defined as a number of SNPs, a
#'   number of base pairs, or a size in centiMorgans (which requires the
#'   'genetic_dist' column of the loci table to be filled in). The statistic is
#'   then recomputed leaving out one block at a time, and the standard error is
#'   computed with the weighted delete-m jackknife of Busing et al. (1999),
#'   where each block is weighted by the number of loci with non-missing values
#'   that it contains. Blocks should be larger than the extent of linkage
#'   disequilibrium (e.g. 5 cM, or a few Mb).
#'
#' @references Busing FMTA, Meijer E, Van Der Leeden R (1999) Delete-m
#'   jackknife for unequal m. Statistics and Computing 9:3-8.
#' @param .x A vector containing the per locus estimates.
#' @param loci_table a dataframe including at least a column 'chromosome', and
#'   additionally a column 'position' if `size_unit` is "bp", or a column
#'   'genetic_dist' if `size_unit` is "cM".
#' @param operator The operator to use to aggregate the per locus estimates.
#'   Either "mean", "sum" or "custom" to use a custom function `f`.
#' @param block_size The size of the blocks to leave out.
#' @param size_unit Either "snp", "cM" or "bp". If "snp", the block size is in
#'   number of SNPs. If "cM", the block size is in centiMorgans. If "bp", the
#'   block size is in base pairs.
#' @param f a custom function to use for the statistic. This function
#'   should take a vector of locus estimates and return a single value. Note
#'   that, for a custom function, the statistic is recomputed from scratch
#'   for every block, which is slower than using "mean" or "sum".
#' @param ... Additional arguments to be passed to the custom operator function.
#' @returns A tibble with columns: 'estimate' (the statistic computed on all
#'   loci), 'jack_estimate' (the bias corrected jackknife estimate), 'se' (the
#'   jackknife standard error), 'z' (the z score, jack_estimate/se), and
#'   'n_blocks' (the number of blocks used).
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' miss_by_locus <- loci_missingness(example_gt)
#'
#' # Calculate the mean missingness and its standard error,
#' # using blocks of 2 SNPs
#' jackknife_stats_generic(miss_by_locus,
#'   loci_table = show_loci(example_gt),
#'   operator = "mean", block_size = 2,
#'   size_unit = "snp"
#' )
#'
jackknife_stats_generic <- function(.x, loci_table,
                                    operator = c("mean", "sum", "custom"),
                                    block_size,
                                    size_unit = c("snp", "cM", "bp"),
                                    f = NULL, ...) {
  size_unit <- match.arg(size_unit)
  operator <- match.arg(operator)
  if (operator == "custom" && is.null(f)) {
    stop("If operator is 'custom', f must be provided.")
  }
  # check the loci_table has same number of rows as x
  if (nrow(loci_table) != length(.x)) {
    stop("loci_table must have the same number of rows as x.")
  }
  # assign each locus to a block
  block_id <- .jackknife_block_ids(
    loci_table = loci_table,
    block_size = block_size,
    size_unit = size_unit
  )

  # get the sums and number of non-missing loci in each block
  x_blocks <- .block_sums(.x, block_id)
  m <- x_blocks$n
  if (operator == "sum") {
    theta_full <- sum(x_blocks$sums)
    theta_loo <- theta_full - x_blocks$sums
  } else if (operator == "mean") {
    theta_full <- sum(x_blocks$sums) / sum(m)
    theta_loo <- (sum(x_blocks$sums) - x_blocks$sums) / (sum(m) - m)
  } else if (operator == "custom") {
    # we have to recompute the statistic leaving out each block in turn
    theta_full <- f(.x, ...)
    theta_loo <- matrix(
      vapply(sort(unique(block_id)), function(i_block) {
        f(.x[block_id != i_block], ...)
      }, numeric(1)),
      ncol = 1
    )
  }
  .jackknife_weighted(
    theta_full = theta_full,
    theta_loo = theta_loo,
    m = m
  )
}
