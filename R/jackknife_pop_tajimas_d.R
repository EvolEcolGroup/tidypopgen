#' Compute Tajima's D with block jackknife standard errors
#'
#' @description This function computes genome-wide Tajima's D, together with
#'   its standard error estimated with a weighted block jackknife.
#'
#' @details The genome is divided into non-overlapping blocks, which never span
#'   more than one chromosome, and Tajima's D is recomputed leaving out one
#'   block at a time. The standard error is computed with the weighted delete-m
#'   jackknife of Busing et al. (1999), where each block is weighted by the
#'   number of loci it contains. Tajima's D is a function of sums across loci
#'   (of the nucleotide diversity, the number of segregating sites, and
#'   quantities that depend on the sample size of each locus), so leaving out
#'   a block is achieved by subtracting its sums, without having to rescan the
#'   genotypes. See [jackknife_stats_generic()] for details.
#'
#'   Missing genotypes are accounted for by using the number of non-missing
#'   allele copies at each locus as its sample size, following Bailey et al.
#'   (2025) as implemented in `pixy`. See [pop_tajimas_d()] for details. Loci
#'   with no genotyped individuals are ignored. The variance of Tajima (1989)
#'   is computed for the mean sample size across all loci (rounded to the
#'   nearest integer), and the same sample size is used when leaving out each
#'   block: recomputing the rounded mean for each block left out would make
#'   Tajima's D jump whenever this changes the rounding (e.g. when the mean
#'   sample size is close to x.5), which the jackknife would amplify into an
#'   inflated standard error and a biased jackknife estimate.
#'
#'   Note that the variance of Tajima (1989) assumes that there is no
#'   recombination among loci. Over a whole genome, Tajima's D is thus not
#'   standardised: it is approximately proportional to
#'   \eqn{\hat{\theta}_\pi / \hat{\theta}_W - 1}, with a factor that depends
#'   on the sample size (about 5 for 10 allele copies, 4.2 for 20, and 3.4
#'   for 100), and with genome-wide data its standard error is much smaller
#'   than 1. Use the jackknife standard error, rather than the rules of thumb
#'   used for single loci, to judge whether genome-wide Tajima's D differs
#'   from 0, and compare groups with different sample sizes with caution. The
#'   z score tests whether genome-wide Tajima's D differs from 0, which most
#'   demographic histories other than a constant population size will cause.
#'   Blocks should be larger than the extent of linkage disequilibrium, or
#'   the standard error will be underestimated.
#'
#' @references Busing FMTA, Meijer E, Van Der Leeden R (1999) Delete-m
#'   jackknife for unequal m. Statistics and Computing 9:3-8.
#'
#'   Tajima F (1989) Statistical method for testing the neutral mutation
#'   hypothesis by DNA polymorphism. Genetics 123:585-595.
#'
#'   Bailey N, Stevison L, Samuk K (2025) Correcting for bias in estimates of
#'   \eqn{\theta_W} and Tajima's D from missing data in next-generation
#'   sequencing. Molecular Ecology Resources, e14104.
#'   \doi{10.1111/1755-0998.14104}
#' @param .x a `gen_tibble` object, optionally grouped by one variable
#' @param type type of object to return, if using grouped method. One of
#'   "matrix", "tidy", or "list". Default is "matrix".
#'   Data can be grouped by more than one variable: groups are then labelled
#'   by the values of the grouping variables pasted with "." (e.g. `pop1.F`),
#'   and the "tidy" output has one column per grouping variable rather than a
#'   single `group` column. The `group` column, or the columns of the
#'   grouping variables, keep the type of the grouping variables (e.g. a
#'   factor keeps its levels), so that groups stay in their order when the
#'   output is sorted or plotted.
#' @param block_size The size of the blocks to leave out.
#' @param size_unit Either "snp", "cM" or "bp". If "snp", the block size is in
#'   number of SNPs. If "cM", the block size is in centiMorgans. If "bp", the
#'   block size is in base pairs.
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns if data is not grouped, a data frame with the following columns:
#' - `estimate`: Tajima's D computed on all loci
#' - `jack_estimate`: the bias corrected jackknife estimate
#' - `se`: the jackknife standard error
#' - `z`: the z score (jack_estimate/se)
#' - `n_blocks`: the number of blocks used
#' if data are grouped, either:
#' a data frame with the following columns:
#' - `quantity`: the quantity reported in the row (`estimate`,
#'   `jack_estimate`, `se`, `z` or `n_blocks`, as described above)
#' - `group`: the value of the quantity for the group (there will be
#'    as many of these columns as groups in the gen_tibble, and they will be
#'    named by the grouping levels)
#' a tidy tibble with the following columns:
#' - `group`: the name of the group
#' - `estimate`, `jack_estimate`, `se`, `z`, `n_blocks`: as described above
#' or a list of data frames, one per group, with the same columns as for
#' ungrouped data.
#' @rdname jackknife_pop_tajimas_d
#' @export
#' @examples
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # Calculate Tajima's D using blocks of 2 SNPs
#' example_gt %>%
#'   jackknife_pop_tajimas_d(
#'     block_size = 2,
#'     size_unit = "snp"
#'   )
#'
jackknife_pop_tajimas_d <- function(.x, ...) {
  UseMethod("jackknife_pop_tajimas_d", .x)
}

#' @export
#' @rdname jackknife_pop_tajimas_d
jackknife_pop_tajimas_d.gen_tbl <- function(
  .x,
  type = c("matrix", "tidy", "list"),
  block_size,
  size_unit = c("snp", "cM", "bp"),
  ...
) {
  rlang::check_dots_empty()
  # Check if the input is a gen_tibble
  stopifnot_gen_tibble(.x)
  type <- match.arg(type)
  size_unit <- match.arg(size_unit)

  # assign each locus to a block
  block_id <- .jackknife_block_ids(
    loci_table = show_loci(.x),
    block_size = block_size,
    size_unit = size_unit
  )
  # sum the quantities needed to compute tajimas d in each block (each block
  # is a non-overlapping window), for each group, in a single pass over the
  # genotypes
  term_sums <- .tajimas_d_sums_by_group(.x,
    locus_first = block_id - 1L,
    locus_last = block_id - 1L,
    n_windows = max(block_id)
  )

  res <- list()
  # now we can loop around each population to compute the jackknife
  for (i_grp in seq_len(ncol(term_sums$k_hat))) {
    res[[i_grp]] <- .jackknife_tajimas_d_blocks(
      .tajimas_d_group_sums(term_sums, i_grp)
    )
  }
  if (length(res) == 1) { # if we only have one pop, return a data.frame
    res <- res[[1]]
    return(res)
  }

  groups <- .group_info(.x)
  names(res) <- groups$labels

  if (type == "matrix") {
    res <- .jackknife_format(dplyr::bind_rows(res),
      stat_names = names(res),
      type = "matrix"
    )
    return(res)
  } else if (type == "tidy") {
    res <- .bind_group_rows(res, groups)
    return(res)
  } else if (type == "list") {
    return(res)
  }
}

#' Weighted block jackknife of Tajima's D from the sums of each block
#'
#' Tajima's D is computed on all loci, and leaving out each block in turn by
#' subtracting its sums. The sample size used in the variance of Tajima (1989)
#' is the rounded mean sample size across all loci, and it is the same for all
#' the leave-one-out estimates: rounding the mean sample size of each
#' leave-one-out estimate would make Tajima's D jump whenever leaving out a
#' block changes the rounding (e.g. when the mean sample size is close to
#' x.5), and the jackknife amplifies these jumps (by about the number of
#' blocks), inflating the standard error and biasing the jackknife estimate.
#'
#' @param term_blocks a matrix with one row per block, and the sums of the
#'   quantities of `tajimas_d_locus_terms()` as columns (`k_hat`, `seg`,
#'   `theta_w`, `n` and `n_loci`), as returned by `.tajimas_d_group_sums()`
#' @returns a tibble as returned by `.jackknife_weighted()`
#' @keywords internal
#' @noRd
.jackknife_tajimas_d_blocks <- function(term_blocks) {
  # sum across all loci, and leaving out each block in turn
  term_tot <- colSums(term_blocks)
  term_loo <- rep_rows(term_tot, nrow(term_blocks)) - term_blocks
  colnames(term_loo) <- colnames(term_blocks)
  # the sample size for the variance, from all loci
  n_var <- tajimas_d_n_var(term_tot[["n"]], term_tot[["n_loci"]])
  # Tajima's D on all loci, and leaving out each block in turn
  .jackknife_weighted(
    theta_full = tajimas_d_from_term_sums(term_tot, n_var = n_var),
    theta_loo = tajimas_d_from_term_sums(term_loo, n_var = n_var),
    # the number of loci with at least one valid allele copy in each block
    m = term_blocks[, "n_loci", drop = FALSE]
  )
}

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
jackknife_pop_tajimas_d.data.frame <- function(.x, ...) {
  jackknife_pop_tajimas_d(
    restore_genotypes_tbl(.x, "jackknife_pop_tajimas_d"),
    ...
  )
}
