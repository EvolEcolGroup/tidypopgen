#' Compute pairwise Fst with block jackknife standard errors
#'
#' @description This function computes genome-wide pairwise Fst, together with
#'   its standard error estimated with a weighted block jackknife.
#'
#' @details The genome is divided into non-overlapping blocks, which never span
#'   more than one chromosome, and Fst is recomputed leaving out one block at a
#'   time. The standard error is computed with the weighted delete-m jackknife
#'   of Busing et al. (1999), where each block is weighted by the number of loci
#'   it contains. Fst is computed as the ratio of the sum of the per locus
#'   numerators and denominators, so leaving out a block is achieved by
#'   subtracting its sums, without having to rescan the genotypes. See
#'   [jackknife_stats_generic()] for details.
#'
#' @references Busing FMTA, Meijer E, Van Der Leeden R (1999) Delete-m
#'   jackknife for unequal m. Statistics and Computing 9:3-8.
#' @param .x a grouped `gen_tibble` object
#' @param type type of object to return. One of "matrix" or "tidy". Default is
#'   "matrix". "matrix" returns a dataframe where each row is a quantity
#'   (estimate, jack_estimate, se, z and n_blocks), followed by columns for
#'   each pairwise population a and b comparison. "tidy" returns a tidy tibble
#'   where each row is one pairwise population a and b comparison.
#' @param method the method to use for calculating Fst, one of 'Hudson',
#'   'Nei87', and 'WC84'. See [pairwise_pop_fst()] for details.
#' @param block_size The size of the blocks to leave out.
#' @param size_unit Either "snp", "cM" or "bp". If "snp", the block size is in
#'   number of SNPs. If "cM", the block size is in centiMorgans. If "bp", the
#'   block size is in base pairs.
#' @returns either a data frame with the following columns:
#' - `quantity`: the quantity reported in the row, one of `estimate` (Fst
#'   computed on all loci), `jack_estimate` (the bias corrected jackknife
#'   estimate), `se` (the jackknife standard error), `z` (jack_estimate/se),
#'   and `n_blocks` (the number of blocks used)
#' - `fst_a.b`: the values for the pairwise Fst between population a and b
#'   (there will be multiple such columns if there are more than two
#'   populations)
#' or a tidy tibble with the following columns:
#' - `stat_name`: the name of population a and b used in the pairwise Fst
#'   calculation (e.g. "fst_pop1.pop2")
#' - `estimate`, `jack_estimate`, `se`, `z`, `n_blocks`: as described above
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>%
#'   group_by(population) %>%
#'   jackknife_pairwise_pop_fst(
#'     block_size = 2,
#'     size_unit = "snp"
#'   )
#'
jackknife_pairwise_pop_fst <- function(.x,
                                       type = c("matrix", "tidy"),
                                       method = c("Hudson", "Nei87", "WC84"),
                                       block_size,
                                       size_unit = c("snp", "cM", "bp")) {
  # Check if the input is a gen_tibble
  stopifnot_gen_tibble(.x)
  method <- match.arg(method)
  type <- match.arg(type)
  size_unit <- match.arg(size_unit)

  # assign each locus to a block
  block_id <- .jackknife_block_ids(
    loci_table = show_loci(.x),
    block_size = block_size,
    size_unit = size_unit
  )
  # compute Fst on all loci, and leaving out each block in turn
  fst_blocks <- .jackknife_fst_blocks(.x, method = method, block_id = block_id)
  # compute the jackknife estimates
  res <- .jackknife_weighted(
    theta_full = fst_blocks$full,
    theta_loo = fst_blocks$loo,
    m = fst_blocks$m
  )
  .jackknife_format(res, stat_names = names(fst_blocks$full), type = type)
}
