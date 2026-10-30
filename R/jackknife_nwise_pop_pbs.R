#' Compute the Population Branch Statistics with block jackknife standard errors
#'
#' @description The function computes the genome-wide population branch
#'   statistics (PBS) for each combination of populations, together with its
#'   standard error estimated with a weighted block jackknife. The PBS is a
#'   measure of the genetic differentiation between one focal population and two
#'   reference populations.
#'
#' @details The genome is divided into non-overlapping blocks, which never span
#'   more than one chromosome, and pairwise Fst (and thus PBS) is recomputed
#'   leaving out one block at a time. The standard error is computed with the
#'   weighted delete-m jackknife of Busing et al. (1999), where each block is
#'   weighted by the number of loci it contains. See [jackknife_stats_generic()]
#'   for details.
#'
#' @references Busing FMTA, Meijer E, Van Der Leeden R (1999) Delete-m
#'   jackknife for unequal m. Statistics and Computing 9:3-8.
#' @param .x a grouped `gen_tibble` object
#' @param type type of object to return. One of "matrix" or "tidy". Default is
#'   "matrix". "matrix" returns a dataframe where each row is a quantity
#'   (estimate, jack_estimate, se, z and n_blocks), followed by columns for
#'   each population comparison. "tidy" returns a tidy tibble where each row is
#'   one population comparison.
#' @param fst_method the method to use for calculating Fst, one of 'Hudson',
#'   'Nei87', and 'WC84'. See [pairwise_pop_fst()] for details.
#' @param return_fst a logical value indicating whether to return the Fst values
#' @param block_size The size of the blocks to leave out.
#' @param size_unit Either "snp", "cM" or "bp". If "snp", the block size is in
#'   number of SNPs. If "cM", the block size is in centiMorgans. If "bp", the
#'   block size is in base pairs.
#' @returns either a data frame with the following columns:
#' - `quantity`: the quantity reported in the row, one of `estimate` (computed
#'   on all loci), `jack_estimate` (the bias corrected jackknife estimate), `se`
#'   (the jackknife standard error), `z` (jack_estimate/se), and `n_blocks` (the
#'   number of blocks used)
#' - `fst_a.b`: the Fst value for population a and b, if `return_fst` is TRUE
#' - `pbs_a.b.c`: the PBS value for population a given b & c (there
#'   will be multiple such columns covering all 3 way combinations of
#'   populations in the grouped `gen_tibble` object)
#' - `pbsn1_a.b.c`: the normalised PBS value for population a given b & c
#' or a tidy tibble with the following columns:
#' - `stat_name`: the name of populations used in the pbs calculation
#'    (e.g. "pbs_pop1.pop2.pop3"). If return_fst is TRUE, stat_name will also
#'    include "fst" calculations in the same column (e.g. "fst_pop1.pop2").
#' - `estimate`, `jack_estimate`, `se`, `z`, `n_blocks`: as described above
#' @export
#' @examples
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # Calculate nwise pbs with blocks of 2 SNPs
#' example_gt %>%
#'   jackknife_nwise_pop_pbs(
#'     block_size = 2,
#'     size_unit = "snp"
#'   )
jackknife_nwise_pop_pbs <- function(.x,
                                    type = c("matrix", "tidy"),
                                    fst_method = c("Hudson", "Nei87", "WC84"),
                                    return_fst = FALSE,
                                    block_size,
                                    size_unit = c("snp", "cM", "bp")) {
  # Check if the input is a grouped gen_tibble
  if (!inherits(.x, "gen_tbl") || !inherits(.x, "grouped_df")) {
    stop(".x should be a grouped gen_tibble")
  }
  # check that we only have one grouping variable
  if (length(.x %>% dplyr::group_vars()) > 1) {
    stop("jackknife_nwise_pop_pbs only works with one grouping variable")
  }
  # Check if the return_fst is logical
  if (!is.logical(return_fst)) {
    stop("return_fst must be a logical value (TRUE or FALSE)")
  }
  type <- match.arg(type)
  size_unit <- match.arg(size_unit)
  # get the populations
  .group_levels <- .x %>% group_keys()
  # Check if there are at least 3 populations
  if (nrow(.group_levels) < 3) {
    stop("At least 3 populations are required to compute PBS.")
  }
  fst_method <- match.arg(fst_method)

  # assign each locus to a block
  block_id <- .jackknife_block_ids(
    loci_table = show_loci(.x),
    block_size = block_size,
    size_unit = size_unit
  )
  # compute Fst on all loci, and leaving out each block in turn
  fst_blocks <- .jackknife_fst_blocks(.x,
    method = fst_method,
    block_id = block_id
  )

  # create all 3 way combinations of populations
  pop_combinations <- utils::combn(.group_levels %>% dplyr::pull(1), 3,
    simplify = FALSE
  )

  # for each combination of populations, compute the pbs on all loci...
  pbs_full <- lapply(pop_combinations,
    pbs_one_triplet,
    fst_values = matrix(fst_blocks$full,
      nrow = 1,
      dimnames = list(NULL, names(fst_blocks$full))
    )
  )
  pbs_full <- unlist(do.call(cbind, pbs_full))
  # ...and leaving out each block in turn
  pbs_loo <- lapply(pop_combinations,
    pbs_one_triplet,
    fst_values = fst_blocks$loo
  )
  pbs_loo <- as.matrix(do.call(cbind, pbs_loo))
  # the number of loci in a block is the minimum across the three pairs
  pbs_m <- lapply(pop_combinations, function(pops) {
    pairs <- paste0("fst_", c(
      paste(pops[1], pops[2], sep = "."),
      paste(pops[1], pops[3], sep = "."),
      paste(pops[2], pops[3], sep = ".")
    ))
    # pbs_one_triplet returns 6 statistics per triplet
    m_triplet <- do.call(
      pmin,
      as.data.frame(fst_blocks$m[, pairs, drop = FALSE])
    )
    matrix(m_triplet, nrow = length(m_triplet), ncol = 6)
  })
  pbs_m <- do.call(cbind, pbs_m)

  # if we return fst, we add it before the pbs values
  if (return_fst) {
    pbs_full <- c(fst_blocks$full, pbs_full)
    pbs_loo <- cbind(fst_blocks$loo, pbs_loo)
    pbs_m <- cbind(fst_blocks$m, pbs_m)
  }
  # compute the jackknife estimates
  res <- .jackknife_weighted(
    theta_full = pbs_full,
    theta_loo = pbs_loo,
    m = pbs_m
  )
  .jackknife_format(res, stat_names = names(pbs_full), type = type)
}
