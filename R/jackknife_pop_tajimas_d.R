#' Compute Tajima's D with block jackknife standard errors
#'
#' @description This function computes genome-wide Tajima's D, together with
#'   its standard error estimated with a weighted block jackknife.
#'
#' @details The genome is divided into non-overlapping blocks, which never span
#'   more than one chromosome, and Tajima's D is recomputed leaving out one
#'   block at a time. The standard error is computed with the weighted delete-m
#'   jackknife of Busing et al. (1999), where each block is weighted by the
#'   number of loci it contains. Tajima's D is a function of the sum of the per
#'   locus nucleotide diversity and of the number of segregating sites, so
#'   leaving out a block is achieved by subtracting its sums, without having to
#'   rescan the genotypes. Loci with missing nucleotide diversity are ignored.
#'   See [jackknife_stats_generic()] for details.
#'
#' @references Busing FMTA, Meijer E, Van Der Leeden R (1999) Delete-m
#'   jackknife for unequal m. Statistics and Computing 9:3-8.
#' @param .x a (potentially grouped) `gen_tibble` object
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
  type <- match.arg(type)
  size_unit <- match.arg(size_unit)

  # assign each locus to a block
  block_id <- .jackknife_block_ids(
    loci_table = show_loci(.x),
    block_size = block_size,
    size_unit = size_unit
  )
  # get the number of sampled allele copies for each group
  n <- .n_allele_copies(.x)
  # get the pi for each locus (if it x is grouped, it will be a list)
  pi_by_locus <- loci_pi(.x, type = "list")

  # recast pi_by_locus as a list of one for just a population
  if (!inherits(pi_by_locus, "list")) {
    pi_by_locus <- list(pi_by_locus)
  }

  res <- list()
  # now we can loop around each population to compute the jackknife
  for (i_grp in seq_along(pi_by_locus)) {
    pi <- pi_by_locus[[i_grp]]
    # sum pi and count segregating sites within each block
    seg <- ifelse(is.na(pi), NA, pi > 0.0 & pi < 1.0)
    pi_blocks <- .block_sums(pi, block_id)
    seg_blocks <- .block_sums(seg, block_id)
    k_tot <- sum(pi_blocks$sums)
    seg_tot <- sum(seg_blocks$sums)
    # Tajima's D on all loci, and leaving out each block in turn
    taj_full <- tajimas_d_from_sums(k_hat = k_tot, seg = seg_tot, n = n[i_grp])
    taj_loo <- tajimas_d_from_sums(
      k_hat = k_tot - pi_blocks$sums,
      seg = seg_tot - seg_blocks$sums,
      n = n[i_grp]
    )
    res[[i_grp]] <- .jackknife_weighted(
      theta_full = taj_full,
      theta_loo = taj_loo,
      m = pi_blocks$n
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
