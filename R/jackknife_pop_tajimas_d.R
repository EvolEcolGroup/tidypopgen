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
#'   with no genotyped individuals are ignored.
#'
#' @references Busing FMTA, Meijer E, Van Der Leeden R (1999) Delete-m
#'   jackknife for unequal m. Statistics and Computing 9:3-8.
#'
#'   Bailey N, Stevison L, Samuk K (2025) Correcting for bias in estimates of
#'   \eqn{\theta_W} and Tajima's D from missing data in next-generation
#'   sequencing. Molecular Ecology Resources, e14104.
#'   \doi{10.1111/1755-0998.14104}
#' @param .x a `gen_tibble` object, optionally grouped by one variable
#' @param type type of object to return, if using grouped method. One of
#'   "matrix", "tidy", or "list". Default is "matrix".
#' @param block_size The size of the blocks to leave out.
#' @param size_unit Either "snp", "cM" or "bp". If "snp", the block size is in
#'   number of SNPs. If "cM", the block size is in centiMorgans. If "bp", the
#'   block size is in base pairs.
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
jackknife_pop_tajimas_d <- function(.x,
                                    type = c("matrix", "tidy", "list"),
                                    block_size,
                                    size_unit = c("snp", "cM", "bp")) {
  # Check if the input is a gen_tibble
  stopifnot_gen_tibble(.x)
  if (length(dplyr::group_vars(.x)) > 1) {
    stop("jackknife_pop_tajimas_d only works with one grouping variable")
  }
  type <- match.arg(type)
  size_unit <- match.arg(size_unit)

  # assign each locus to a block
  block_id <- .jackknife_block_ids(
    loci_table = show_loci(.x),
    block_size = block_size,
    size_unit = size_unit
  )
  # get pi and the number of valid allele copies for each locus (one column
  # per group)
  pi_n <- .pi_n_by_group(.x)

  res <- list()
  # now we can loop around each population to compute the jackknife
  for (i_grp in seq_len(ncol(pi_n$pi))) {
    # get the quantities needed to compute tajimas d at each locus
    locus_terms <- tajimas_d_locus_terms(pi_n$pi[, i_grp], pi_n$n[, i_grp])
    # sum each of these quantities within each block
    term_blocks <- .block_sums(locus_terms, block_id)
    # sum across all loci, and leaving out each block in turn
    term_tot <- colSums(term_blocks$sums)
    term_loo <- rep_rows(term_tot, nrow(term_blocks$sums)) - term_blocks$sums
    colnames(term_loo) <- colnames(locus_terms)
    # Tajima's D on all loci, and leaving out each block in turn
    taj_full <- tajimas_d_from_term_sums(term_tot)
    taj_loo <- tajimas_d_from_term_sums(term_loo)
    res[[i_grp]] <- .jackknife_weighted(
      theta_full = taj_full,
      theta_loo = taj_loo,
      # the number of loci with valid pi in each block
      m = term_blocks$n[, 1, drop = FALSE]
    )
  }
  if (length(res) == 1) { # if we only have one pop, return a data.frame
    res <- res[[1]]
    return(res)
  }

  names(res) <- dplyr::group_keys(.x) %>% pull(1)

  if (type == "matrix") {
    res <- .jackknife_format(dplyr::bind_rows(res),
      stat_names = names(res),
      type = "matrix"
    )
    return(res)
  } else if (type == "tidy") {
    res <- bind_rows(res, .id = "group")
    return(res)
  } else if (type == "list") {
    return(res)
  }
}
