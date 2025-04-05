#' Compute pairwise Fst for a sliding window
#'
#' @description
#' This function computes pairwise Fst for a sliding window across a set of loci. It uses the Hudson method for calculating Fst.
#'
#' @param x a `gen_tibble` object
#' @param window_size the size of the sliding window (in number of loci)
#' @param step_size the step size for the sliding window (in number of loci)
#' @param min_loci the minimum number of loci required in a window to compute Fst
#' @returns a data frame with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `window_start`: the starting locus of the window
#' - `window_end`: the ending locus of the window
#' - `fst_a_b`: the pairwise Fst value for the population a and b (there will be
#' multiple such columns if there are more than two populations)
#' @export

window_pairwise_pop_fst <- function(x, window_size = 10, step_size = 1, min_loci = 5) {
  # Check if the input is a gen_tibble
  stopifnot_gen_tibble(x)

  # create the pairwise Fst by locus, saving numerator and denominator
  pair_fst <- pairwise_pop_fst(x, method = "Hudson", return_num_dem = TRUE)

  # Get the loci information
  loci_info <- show_loci(x)
}
