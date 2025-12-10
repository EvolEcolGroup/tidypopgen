#' Compute Tajima's D for a sliding window
#'
#' @description This function computes Tajima's D for a sliding window across
#'   each chromosome.
#'
#' @param .x a (potentially grouped) `gen_tibble` object
#' @param type type of object to return, if using grouped method. One of
#'   "matrix", "tidy", or "list". Default is "matrix".
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
#' @returns if data is not grouped, a data frame with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `start`: the starting locus of the window
#' - `end`: the ending locus of the window
#' - `tajimas_d`: the Tajima's D for the population
#' if data are grouped, either:
#' a data frame as above with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `start`: the starting locus of the window
#' - `end`: the ending locus of the window
#' - `n_loci`: the number of loci in the window
#' - `group`: the Tajima's D for the group for the given window  (there will be
#'    as many of these columns as groups in the gen_tibble, and they will be
#'    named by the grouping levels)
#' a tidy tibble with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `start`: the starting locus of the window
#' - `end`: the ending locus of the window
#' - `n_loci`: the number of loci in the window
#' - `group`: the name of the group
#' - `stat`: the Tajima's D for the given group at the given window
#' or a list of data frames, one per group, with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `start`: the starting locus of the window
#' - `end`: the ending locus of the window
#' - `stat`: the Tajima's D for the given window
#' - `n_loci`: the number of loci in the window
#' @export
#' @examples
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # Calculate Tajima's D across a window of 3 SNPs, with a step size of 2 SNPs
#' example_gt %>%
#'   windows_pop_tajimas_d(
#'     window_size = 3, step_size = 2,
#'     size_unit = "snp", min_loci = 2
#'   )
#'
windows_pop_tajimas_d <- function(.x,
                                  type = c("matrix", "tidy", "list"),
                                  window_size,
                                  step_size,
                                  size_unit = c("snp", "bp"),
                                  min_loci = 1,
                                  complete = FALSE) {
  # Check if the input is a gen_tibble
  stopifnot_gen_tibble(.x)
  type <- match.arg(type)

  # if x is grouped, get the pop sizes for each group
  if (inherits(.x, "grouped_gen_tbl")) {
    # get the population sizes
    n <- .x %>%
      dplyr::summarise(n = n()) %>%
      dplyr::pull(.data$n)
  } else {
    # if not grouped, just use the number of individuals
    n <- nrow(.x)
  }
  # get the pi for each locus (if it x is grouped, it will be a list)
  pi_by_locus <- loci_pi(.x, type = "list")

  # recast pi_by_locus as a list of one for just a population
  if (!inherits(pi_by_locus, "list")) {
    pi_by_locus <- list(pi_by_locus)
  }

  res <- list()
  # now we can loop around each population to compute the windows
  for (i_grp in seq_along(pi_by_locus)) {
    window_taj <- windows_stats_generic(
      .x = pi_by_locus[[i_grp]],
      loci_table = show_loci(.x),
      operator = "custom",
      window_size = window_size,
      step_size = step_size,
      size_unit = size_unit,
      min_loci = min_loci,
      complete = complete,
      f = tajimas_d_from_pi_vec,
      n = n[i_grp] * 2 # because we need the number of alleles
    )
    res[[i_grp]] <- window_taj
  }
  if (length(res) == 1) { # if we only have one pop, return a data.frame
    res <- res[[1]]
    return(res)
  }

  names(res) <- dplyr::group_keys(.x) %>% pull(1)

  if (type == "matrix") {
    res <- bind_rows(res, .id = "group")
    res <-
      res %>% tidyr::pivot_wider(names_from = "group", values_from = "stat")
    return(res)
  } else if (type == "tidy") {
    res <- bind_rows(res, .id = "group")
    return(res)
  } else if (type == "list") {
    return(res)
  }
}
