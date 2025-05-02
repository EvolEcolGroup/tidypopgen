#' Compute the Population Branch Statistics for each combination of populations
#'
#' The function computes the population branch statistics (PBS) for each
#' combination of populations at each locus. The PBS is a measure of the genetic
#' differentiation between one focal population and two reference populations,
#' and is used to identify outlier loci that may be under selection.
#' @references Yi X, et al. (2010) Sequencing of 50 human exomes reveals
#' adaptation to high altitude. Science 329: 75-78.
#' @param .x A grouped `gen_tibble`
#' @param fst_method A character string specifying the method to use for
#'   computing Fst
#' @param return_fst A logical value indicating whether to return the Fst values
#'   along with the PBS values. Default is `FALSE`.
#' @return A tibble with the following columns:
#' - `pbs_a.b.c`: the PBS value for population a given b & c (there
#'   will be multiple such columns covering all 3 way combinations of
#'   populations in the grouped `gen_tibble` object)
#' @export

pairwise_pop_pbs <- function(.x, fst_method = c("Hudson"), return_fst = FALSE) {
  # Check if the input is a grouped gen_tibble
  if (!inherits(.x, "gen_tbl") || !inherits(.x, "grouped_df")) {
    stop(".x should be a grouped gen_tibble")
  }
  # @TODO check that we have only one grouping variable.

  # get the populations
  .group_levels <- .x %>% group_keys()
  # Check if there are at least 3 populations
  if (nrow(.group_levels) < 3) {
    stop("At least 3 populations are required to compute PBS.")
  }

  # Compute pairwise Fst values
  fst_values <- pairwise_pop_fst(.x,
    method = fst_method,
    by_locus = TRUE
  )$Fst_by_locus

  # create all 3 way combinations of populations
  pop_combinations <- utils::combn(.group_levels %>% dplyr::pull(1), 3,
    simplify = FALSE
  )

  # for each combination of populations, compute the pbs
  pbs_results <- lapply(pop_combinations,
    pbs_one_triplet,
    fst_values = fst_values
  )
  # combine all tables
  pbs_results <- do.call(cbind, pbs_results)
  # set row names same as the fst table
  rownames(pbs_results) <- rownames(fst_values)
  # add fst values if requested
  if (return_fst) {
    pbs_results <- cbind(fst_values, pbs_results)
  }
  return(pbs_results)
}


#' Internal function to compute pbs for a single triplet of populations
#'
#' @param pops A vector of 3 population names
#' @param fst_values A matrix of Fst values (rows are markers/windows, columns
#' are named as fst_pop1.pop2, fst_pop1.pop3, fst_pop2.pop3)
#' @return A tibble with the PBS values for the triplet of populations
#' @keywords internal

pbs_one_triplet <- function(pops, fst_values) {
  pop1 <- pops[1]
  pop2 <- pops[2]
  pop3 <- pops[3]

  # Extract Fst values for the current combination as vectors
  # if we have a tibble (as obtained from windows_pairwise_pop_fst), we need to
  # pull the columns with the fst values
  if (inherits(fst_values, "tbl_df")) {
    fst12 <- fst_values %>% dplyr::pull(paste0("fst_", pop1, ".", pop2))
    fst13 <- fst_values %>% dplyr::pull(paste0("fst_", pop1, ".", pop3))
    fst23 <- fst_values %>% dplyr::pull(paste0("fst_", pop2, ".", pop3))
  } else { # if we have a matrix, we just subset the columns
    fst12 <- fst_values[, paste0("fst_", pop1, ".", pop2)]
    fst13 <- fst_values[, paste0("fst_", pop1, ".", pop3)]
    fst23 <- fst_values[, paste0("fst_", pop2, ".", pop3)]
  }

  # Compute branch lengths
  t12 <- -log(1 - fst12)
  t13 <- -log(1 - fst13)
  t23 <- -log(1 - fst23)

  # Compute PBS values
  pbs_1 <- (t12 + t13 - t23) / 2
  pbs_2 <- (t12 + t23 - t13) / 2
  pbs_3 <- (t13 + t23 - t12) / 2

  pbsn1_1 <- pbs_1 / (1 + pbs_1 + pbs_2 + pbs_3)
  pbsn1_2 <- pbs_2 / (1 + pbs_1 + pbs_2 + pbs_3)
  pbsn1_3 <- pbs_3 / (1 + pbs_1 + pbs_2 + pbs_3)

  # Create a tibble with the results
  tibble(
    !!paste0("pbs_", paste(pop1, pop2, pop3, sep = ".")) := pbs_1,
    !!paste0("pbs_", paste(pop2, pop1, pop3, sep = ".")) := pbs_2,
    !!paste0("pbs_", paste(pop3, pop1, pop2, sep = ".")) := pbs_3,
    !!paste0("pbsn1_", paste(pop1, pop2, pop3, sep = ".")) := pbsn1_1,
    !!paste0("pbsn1_", paste(pop2, pop1, pop3, sep = ".")) := pbsn1_2,
    !!paste0("pbsn1_", paste(pop3, pop1, pop2, sep = ".")) := pbsn1_3,
  )
}
