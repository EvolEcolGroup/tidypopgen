#' Draw the blocks of bootstrap replicates
#'
#' This function is designed for developers of packages that add bootstrap
#' methods to the generics of tidypopgen (e.g. `tidypopgl`, for genotype
#' likelihoods), so that their replicates resample the same blocks as the
#' methods of tidypopgen. The loci are divided into non-overlapping blocks
#' that never span chromosomes, defined as in [jackknife_stats_generic()] (a
#' number of SNPs, of base pairs, or a size in centiMorgans), and each
#' replicate draws as many blocks as there are, with replacement.
#'
#' The blocks are drawn with the random number generator of R: set its seed
#' (e.g. with [set.seed()]) for reproducible replicates. The matrix of weights
#' has one value per replicate and block; to bound its memory, its size is
#' capped at 1e9 values (option `tidypopgen.sfs_max_bootstrap_weights`).
#'
#' @param loci_table a dataframe including at least a column 'chromosome', and
#'   additionally a column 'position' if `size_unit` is "bp", or a column
#'   'genetic_dist' if `size_unit` is "cM", with one row per locus, in the
#'   order in which the loci are used.
#' @param n_bootstrap the number of replicates.
#' @param block_size the size of the blocks, in the unit of `size_unit`.
#' @param size_unit either "snp", "cM" or "bp". If "snp", the block size is
#'   in number of loci. If "cM", the block size is in centiMorgans. If "bp",
#'   the block size is in base pairs.
#' @returns a list with `block_id` (the block of each locus, an integer from
#'   1 to the number of blocks, numbered in their order of appearance) and
#'   `weights` (an integer matrix with the number of times each block, in
#'   columns, is drawn in each replicate, in rows).
#' @keywords internal
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' set.seed(1)
#' blocks <- bootstrap_blocks(show_loci(example_gt),
#'   n_bootstrap = 3,
#'   block_size = 2, size_unit = "snp"
#' )
#' blocks$block_id
#' blocks$weights
bootstrap_blocks <- function(loci_table, n_bootstrap, block_size,
                             size_unit = c("snp", "cM", "bp")) {
  size_unit <- match.arg(size_unit)
  .sfs_check_bootstrap_args(n_bootstrap, block_size)
  # the blocks of the block jackknife
  block_id <- .jackknife_block_ids(
    loci_table = loci_table,
    block_size = block_size,
    size_unit = size_unit
  )
  n_blocks <- max(block_id)
  if (n_blocks < 2) {
    stop(
      "the block bootstrap needs at least 2 blocks; use a smaller ",
      "block_size"
    )
  }
  # the weights are a matrix with one value per replicate and block, which
  # is copied (as doubles) to project the replicates: with single SNPs as
  # blocks (block_size = 1) of tens of millions of loci, it would need tens
  # of GB, so we cap its size (with an option, for testing)
  max_weights <- getOption("tidypopgen.sfs_max_bootstrap_weights", 1e9)
  n_weights <- as.numeric(n_blocks) * n_bootstrap
  if (n_weights > max_weights) {
    stop(
      "the bootstrap would need ", format(n_weights), " weights (",
      n_bootstrap, " replicates of ", n_blocks, " blocks), which would use ",
      "too much memory; use a larger block_size or fewer replicates"
    )
  }
  weights <- t(stats::rmultinom(n_bootstrap,
    size = n_blocks,
    prob = rep(1, n_blocks)
  ))
  list(block_id = block_id, weights = weights)
}
