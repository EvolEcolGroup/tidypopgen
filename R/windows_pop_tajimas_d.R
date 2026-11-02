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

  # Get the number of sampled allele copies for each group. Pseudohaploid
  # genotypes are stored as 0/2 dosages but contribute one allele per sample.
  if (inherits(.x, "grouped_gen_tbl")) {
    n <- .x %>%
      dplyr::summarise(n = sum(indiv_ploidy(.data$genotypes))) %>%
      dplyr::pull(.data$n)
  } else {
    n <- sum(indiv_ploidy(.x))
  }
  # get the pi for each locus (if it x is grouped, it will be a list)
  pi_by_locus <- loci_pi(.x, type = "list")

  # recast pi_by_locus as a list of one for just a population
  if (!inherits(pi_by_locus, "list")) {
    pi_by_locus <- list(pi_by_locus)
  }

  # define the windows once (as in windows_stats_generic())
  size_unit <- match.arg(size_unit)
  loci_table <- show_loci(.x)
  .check_window_args(
    loci_table,
    n_values = nrow(loci_table),
    window_size = window_size,
    step_size = step_size,
    size_unit = size_unit,
    min_loci = min_loci,
    complete = complete
  )
  win <- .window_definitions(
    loci_table,
    window_size = window_size,
    step_size = step_size,
    size_unit = size_unit
  )
  # windows without loci are NA
  empty <- win$hi == win$lo

  res <- list()
  # now we can loop around each population to compute the windows
  for (i_grp in seq_along(pi_by_locus)) {
    pi <- as.numeric(pi_by_locus[[i_grp]])
    # tajimas d only needs, for each window, the sum of pi (accumulated in
    # long double, as sum() does) and the number of segregating sites
    win_sums <- window_sums_cpp(
      matrix(c(pi, as.numeric(pi > 0)), ncol = 2),
      loci = win$loci,
      lo = win$lo,
      hi = win$hi,
      long_double = TRUE
    )
    stat <- tajimas_d_from_k_hat_seg(
      k_hat = win_sums$sum[, 1],
      seg = win_sums$sum[, 2],
      n = n[i_grp]
    )
    # the number of loci with a non-missing pi in each window
    n_loci <- win_sums$n[, 1]
    stat[empty] <- NA
    n_loci[empty] <- NA
    # if we only want complete windows, set incomplete ones to NA
    if (complete) {
      stat[!win$complete] <- NA
      n_loci[!win$complete] <- NA
    }
    window_taj <- win$windows
    window_taj$stat <- stat
    window_taj$n_loci <- n_loci
    # set values to NA if there are fewer than min_loci loci
    window_taj$stat[window_taj$n_loci < min_loci] <- NA
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
