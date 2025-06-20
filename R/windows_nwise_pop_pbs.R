#' Compute the Population Branch Statistics over a sliding window
#'
#' @description The function computes the population branch statistics (PBS) for
#'   a sliding window for each combination of populations at each locus. The
#'   PBS is a measure of the genetic differentiation between one focal
#'   population and two reference populations, and is used to identify outlier
#'   loci that may be under selection.
#'
#' @param .x a grouped `gen_tibble` object
#' @param type type of object to return. One of "matrix" or "tidy". Default is
#'   "matrix". "matrix" returns a dataframe where each row is a window, followed
#'   by columns of pbs values for each population comparison. "tidy" returns a
#'   tidy tibble of the same data in 'long' format, where each
#'   row is one window for one population comparison.
#' @param fst_method the method to use for calculating Fst, one of 'Hudson',
#'   'Nei87', and 'WC84'. See [pairwise_pop_fst()] for details.
#' @param return_fst a logical value indicating whether to return the Fst values
#' @param window_size The size of the window to use for the estimates.
#' @param step_size The step size to use for the windows.
#' @param size_unit Either "snp" or "bp". If "snp", the window size and step
#'   size are in number of SNPs. If "bp", the window size and step size are in
#'   base pairs.
#' @param min_loci The minimum number of loci required to calculate a window
#'   statistic. If the number of loci in a window is less than this, the window
#'   statistic will be NA.
#' @param complete Should the function be evaluated on complete windows only? If
#'   FALSE, the default, then partial computations will be allowed at the end of
#'   the chromosome.
#' @returns either a data frame with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `start`: the starting locus of the window
#' - `end`: the ending locus of the window
#' - `pbs_a.b.c`: the PBS value for population a given b & c (there
#'   will be multiple such columns covering all 3 way combinations of
#'   populations in the grouped `gen_tibble` object)
#' - `fst_a.b`: the Fst value for population a and b, if `return_fst` is TRUE
#' or a tidy tibble with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `start`: the starting locus of the window
#' - `end`: the ending locus of the window
#' - `stat_name`: the name of populations used in the pbs calculation
#'    (e.g. "pbs_pop1.pop2.pop3"). If return_fst is TRUE, stat_name will also
#'    include "fst" calculations in the same column (e.g. "fst_pop1.pop2").
#' - `value`: the pbs value for the populations
#' @export
#' @examples
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # Calculate nwise pbs across a window of 3 SNPs, with a step size of 2 SNPs
#' example_gt %>%
#'   windows_nwise_pop_pbs(
#'     window_size = 3, step_size = 2,
#'     size_unit = "snp", min_loci = 2
#'   )
windows_nwise_pop_pbs <- function(.x,
                                  type = c("matrix", "tidy"),
                                  fst_method = c("Hudson", "Nei87", "WC84"),
                                  return_fst = FALSE,
                                  window_size,
                                  step_size,
                                  size_unit = c("snp", "bp"),
                                  min_loci = 1,
                                  complete = FALSE) {
  # Check if the input is a grouped gen_tibble
  if (!inherits(.x, "gen_tbl") || !inherits(.x, "grouped_df")) {
    stop(".x should be a grouped gen_tibble")
  }
  type <- match.arg(type)
  # get the populations
  .group_levels <- .x %>% group_keys()
  # Check if there are at least 3 populations
  if (nrow(.group_levels) < 3) {
    stop("At least 3 populations are required to compute PBS.")
  }
  fst_method <- match.arg(fst_method)

  # create the pairwise Fst by locus, saving numerator and denominator
  fst_values <- windows_pairwise_pop_fst(.x,
    method = fst_method,
    window_size = window_size,
    step_size = step_size,
    size_unit = size_unit,
    min_loci = min_loci,
    complete = complete
  )

  # check that we only have one grouping variable
  if (length(.x %>% dplyr::group_vars()) > 1) {
    stop("windows_nwise_pop_pbs only works with one grouping variable")
  }

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
  # if we returns fst, we simply cbind the whole thing
  if (return_fst) {
    pbs_results <- dplyr::bind_cols(fst_values, pbs_results)
  } else {
    # we just bind the window information
    pbs_results <- dplyr::bind_cols(
      fst_values %>%
        dplyr::select(dplyr::all_of(c("chromosome", "start", "end"))),
      pbs_results
    )
  }

  if (type == "matrix") {
    return(pbs_results)
  } else if (type == "tidy") {
    cols <-
      names(pbs_results)[names(pbs_results) != c("chromosome", "start", "end")]
    pbs_results <-
      pbs_results %>%
      tidyr::pivot_longer(cols = all_of(cols), names_to = "stat_name")
    return(pbs_results)
  }
}
