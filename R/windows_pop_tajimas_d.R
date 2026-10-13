#' Compute Tajima's D for a sliding window
#'
#' @description This function computes Tajima's D for a sliding window across
#'   each chromosome.
#'
#' @details Missing genotypes are accounted for by using the number of
#'   non-missing allele copies at each locus as its sample size, following
#'   Bailey et al. (2025) as implemented in `pixy`. See [pop_tajimas_d()] for
#'   details. Loci with no genotyped individuals are ignored, and are not
#'   counted in `n_loci`. With windows defined in base pairs and
#'   `step_size` equal to `window_size`, the windows match those used by
#'   `pixy`.
#'
#' @references Bailey N, Stevison L, Samuk K (2025) Correcting for bias in
#'   estimates of \eqn{\theta_W} and Tajima's D from missing data in
#'   next-generation sequencing. Molecular Ecology Resources, e14104.
#'   \doi{10.1111/1755-0998.14104}
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

  # get pi and the number of valid allele copies for each locus (one column
  # per group)
  pi_n <- .pi_n_by_group(.x)

  res <- list()
  # now we can loop around each population to compute the windows
  for (i_grp in seq_len(ncol(pi_n$pi))) {
    # get the quantities needed to compute tajimas d at each locus
    locus_terms <- tajimas_d_locus_terms(pi_n$pi[, i_grp], pi_n$n[, i_grp])
    # sum each of these quantities within each window
    window_sums <- list()
    for (i_term in colnames(locus_terms)) {
      window_sums[[i_term]] <- windows_stats_generic(
        .x = locus_terms[, i_term],
        loci_table = show_loci(.x),
        operator = "sum",
        window_size = window_size,
        step_size = step_size,
        size_unit = size_unit,
        min_loci = min_loci,
        complete = complete
      )
    }
    # use the windows for pi as the template for the results (n_loci in this
    # table is the number of loci with at least one valid allele copy in each
    # window)
    window_taj <- window_sums$k_hat
    # compute tajimas d for each window from the sums
    window_taj$stat <- tajimas_d_from_sums(
      k_hat = window_sums$k_hat$stat,
      seg = window_sums$seg$stat,
      theta_w = window_sums$theta_w$stat,
      n_sum = window_sums$n$stat,
      n_loci = window_sums$n_loci$stat
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
