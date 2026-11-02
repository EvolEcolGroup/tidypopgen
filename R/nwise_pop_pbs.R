#' Compute the Population Branch Statistics for each combination of populations
#'
#' The function computes the population branch statistics (PBS) for each
#' combination of populations at each locus. The PBS is a measure of the genetic
#' differentiation between one focal population and two reference populations,
#' and is used to identify outlier loci that may be under selection.
#' @references Yi X, et al. (2010) Sequencing of 50 human exomes reveals
#'   adaptation to high altitude. Science 329: 75-78.
#' @param .x A grouped `gen_tibble`
#' @param type type of object to return. One of "tidy" or "matrix".
#'   Default is "tidy".
#' @param fst_method the method to use for calculating Fst, one of 'Hudson',
#'   'Nei87', and 'WC84'. See [pairwise_pop_fst()] for details.
#' @param return_fst A logical value indicating whether to return the Fst values
#'   along with the PBS values. Default is `FALSE`.
#' @return Either a [matrix_int_names] with locus ID as row names (use
#'   [row_names()] to retrieve them, as they can be integers) and the following
#'   columns:
#' - `pbs_a.b.c`: the PBS value for population a given b & c (there
#'   will be multiple such columns covering all 3 way combinations of
#'   populations in the grouped `gen_tibble` object)
#' - `pbsn1_a.b.c`: the normalized PBS value for population a given b & c.
#' - `fst_a.b`: the Fst value for population a and b, if `return_fst` is TRUE
#' or a tidy tibble with the following columns:
#'  - `loci`: the locus ID
#'  - `stat_name`: the name of populations used in the pbs calculation
#'    (e.g. "pbs_pop1.pop2.pop3"). If return_fst is TRUE, stat_name will also
#'    include "fst" calculations in the same column (e.g. "fst_pop1.pop2").
#' - `value`: the pbs value for the populations
#' @export
#' @examples
#' example_gt <- load_example_gt()
#'
#' # We can compute the PBS for all populations using "Hudson" method
#' example_gt %>%
#'   group_by(population) %>%
#'   nwise_pop_pbs(fst_method = "Hudson")
nwise_pop_pbs <- function(.x,
                          type = c("tidy", "matrix"),
                          fst_method = c("Hudson", "Nei87", "WC84"),
                          return_fst = FALSE) {
  # Check if the input is a grouped gen_tibble
  if (!inherits(.x, "gen_tbl") || !inherits(.x, "grouped_df")) {
    stop(".x should be a grouped gen_tibble")
  }
  # check that we only have one grouping variable
  if (length(.x %>% dplyr::group_vars()) > 1) {
    stop("nwise_pop_pbs only works with one grouping variable")
  }
  # Check if the return_fst is logical
  if (!is.logical(return_fst)) {
    stop("return_fst must be a logical value (TRUE or FALSE)")
  }
  type <- match.arg(type)
  fst_method <- match.arg(fst_method)

  # get the populations
  .group_levels <- .x %>% group_keys()
  # Check if there are at least 3 populations
  if (nrow(.group_levels) < 3) {
    stop("At least 3 populations are required to compute PBS.")
  }

  # Compute pairwise Fst values
  fst_values <- pairwise_pop_fst(.x,
    method = fst_method,
    by_locus = TRUE,
    by_locus_type = "matrix"
  )$Fst_by_locus

  # get the loci names (integer or character) from the fst matrix
  loci_ids <- row_names(fst_values)
  # the fst values as a plain matrix (loci as rows, pairs as columns)
  fst_values <- matrix(
    as.numeric(fst_values),
    nrow = nrow(fst_values),
    dimnames = list(NULL, col_names(fst_values))
  )
  # compute the pbs for all 3 way combinations of populations, in a
  # single matrix (loci as rows)
  pbs_results <- pbs_all_triplets(
    fst_values,
    pops = .group_levels %>% dplyr::pull(1)
  )
  # add fst values if requested
  if (return_fst) {
    pbs_results <- cbind(fst_values, pbs_results)
  }

  if (type == "matrix") {
    pbs_results <- matrix_int_names(pbs_results, row_names = loci_ids)
    return(pbs_results)
  } else if (type == "tidy") {
    # build the long tibble directly (one row per locus and statistic, with
    # the statistics of each locus in consecutive rows), rather than
    # pivoting a data.frame
    pbs_results <- tibble::tibble(
      loci = rep(loci_ids, each = ncol(pbs_results)),
      stat_name = rep(colnames(pbs_results), times = nrow(pbs_results)),
      value = as.vector(t(pbs_results))
    )
    return(pbs_results)
  }
}


#' Internal function to compute pbs for all triplets of populations
#'
#' The branch length \eqn{T = -log(1 - Fst)} of each pair of populations is
#' computed once, and the PBS (and normalised PBS) of each triplet are written
#' into a preallocated matrix.
#'
#' @param fst_values A numeric matrix of Fst values (rows are markers/windows,
#'   columns are named as fst_pop1.pop2, fst_pop1.pop3, fst_pop2.pop3, etc.)
#' @param pops A vector of population names (in the order used to create the
#'   pairs)
#' @return A matrix with one row per row of `fst_values`, and 6 columns for
#'   each triplet of populations a, b and c (as given by `utils::combn()`):
#'   `pbs_a.b.c`, `pbs_b.a.c`, `pbs_c.a.b`, `pbsn1_a.b.c`, `pbsn1_b.a.c` and
#'   `pbsn1_c.a.b`
#' @keywords internal
#' @noRd
pbs_all_triplets <- function(fst_values, pops) {
  # branch lengths for all pairs of populations (computed once per pair)
  t_values <- -log(1 - fst_values)
  # create all 3 way combinations of populations
  pop_combinations <- utils::combn(pops, 3, simplify = FALSE)
  # preallocate the results
  pbs_names <- unlist(lapply(pop_combinations, function(pops) {
    c(
      paste0("pbs_", paste(pops[1], pops[2], pops[3], sep = ".")),
      paste0("pbs_", paste(pops[2], pops[1], pops[3], sep = ".")),
      paste0("pbs_", paste(pops[3], pops[1], pops[2], sep = ".")),
      paste0("pbsn1_", paste(pops[1], pops[2], pops[3], sep = ".")),
      paste0("pbsn1_", paste(pops[2], pops[1], pops[3], sep = ".")),
      paste0("pbsn1_", paste(pops[3], pops[1], pops[2], sep = "."))
    )
  }))
  pbs_results <- matrix(NA_real_,
    nrow = nrow(fst_values),
    ncol = length(pbs_names),
    dimnames = list(NULL, pbs_names)
  )
  for (i_trip in seq_along(pop_combinations)) {
    pop1 <- pop_combinations[[i_trip]][1]
    pop2 <- pop_combinations[[i_trip]][2]
    pop3 <- pop_combinations[[i_trip]][3]
    # branch lengths for the three pairs of this triplet
    t12 <- t_values[, paste0("fst_", pop1, ".", pop2)]
    t13 <- t_values[, paste0("fst_", pop1, ".", pop3)]
    t23 <- t_values[, paste0("fst_", pop2, ".", pop3)]

    # Compute PBS values
    pbs_1 <- (t12 + t13 - t23) / 2
    pbs_2 <- (t12 + t23 - t13) / 2
    pbs_3 <- (t13 + t23 - t12) / 2
    # normalisation term
    pbs_sum <- 1 + pbs_1 + pbs_2 + pbs_3

    # fill the 6 columns of this triplet
    first_col <- 6 * (i_trip - 1)
    pbs_results[, first_col + 1] <- pbs_1
    pbs_results[, first_col + 2] <- pbs_2
    pbs_results[, first_col + 3] <- pbs_3
    pbs_results[, first_col + 4] <- pbs_1 / pbs_sum
    pbs_results[, first_col + 5] <- pbs_2 / pbs_sum
    pbs_results[, first_col + 6] <- pbs_3 / pbs_sum
  }
  pbs_results
}
