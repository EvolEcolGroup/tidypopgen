#' Compute the Population Branch Statistics over a sliding window
#'
#' @description The function computes the population branch statistics (PBS) for
#'   a sligiding window for each combination of populations at each locus. The
#'   PBS is a measure of the genetic differentiation between one focal
#'   population and two reference populations, and is used to identify outlier
#'   loci that may be under selection.
#'
#' @param x a grouped `gen_tibble` object
#' @param method the method to use for calculating Fst. Currently only "Hudson"
#'   is supported.
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
#' @returns a data frame with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `start`: the starting locus of the window
#' - `end`: the ending locus of the window
#' - `pbs_a.b.c`: the PBS value for population a given b & c (there
#'   will be multiple such columns covering all 3 way combinations of
#'   populations in the grouped `gen_tibble` object)
#' - `fst_a.b`: the Fst value for population a and b, if `return_fst` is TRUE
#' @export

window_pairwise_pop_pbs <- function(.x,
                                    method = "Hudson",
                                    return_fst = FALSE,
                                    window_size,
                                    step_size,
                                    size_unit = c("snp", "bp"),
                                    min_loci = 1,
                                    complete = FALSE) {
  message("This is a new function and not fully tested; use it with care")
  # Check if the input is a gen_tibble
  # Check if the input is a grouped gen_tibble
  if (!inherits(.x, "gen_tbl") || !inherits(.x, "grouped_df")) {
    stop(".x should be a grouped gen_tibble")
  }
  # get the populations
  .group_levels <- .x %>% group_keys()
  # Check if there are at least 3 populations
  if (nrow(.group_levels) < 3) {
    stop("At least 3 populations are required to compute PBS.")
  }
  method <- match.arg(method)

  # create the pairwise Fst by locus, saving numerator and denominator
  fst_values <- window_pairwise_pop_fst(.x,
    method = method,
    window_size = window_size,
    step_size = step_size,
    size_unit = size_unit,
    min_loci = min_loci,
    complete = complete
  )

  # create all 3 way combinations of populations
  pop_combinations <- combn(.group_levels %>% dplyr::pull(1), 3,
                            simplify = FALSE)
  # @TODO the above will only work if we have one grouping level

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
  return(pbs_results)
}
