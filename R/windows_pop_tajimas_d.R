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
#'   Data can be grouped by more than one variable: groups are then labelled
#'   by the values of the grouping variables pasted with "." (e.g. `pop1.F`),
#'   and the "tidy" output has one column per grouping variable rather than a
#'   single `group` column. The `group` column, or the columns of the
#'   grouping variables, keep the type of the grouping variables (e.g. a
#'   factor keeps its levels), so that groups stay in their order when the
#'   output is sorted or plotted.
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
#' - `n_loci`: the number of loci in the window (including loci with no
#'    genotyped individuals in some groups; the number of loci with data in
#'    each group is given by the "tidy" and "list" outputs)
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
#' @param ... other arguments passed to specific methods, currently
#'   unused.
#' @rdname windows_pop_tajimas_d
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
windows_pop_tajimas_d <- function(.x, ...) {
  UseMethod("windows_pop_tajimas_d", .x)
}

#' @export
#' @rdname windows_pop_tajimas_d
windows_pop_tajimas_d.gen_tbl <- function(
  .x,
  type = c("matrix", "tidy", "list"),
  window_size,
  step_size,
  size_unit = c("snp", "bp"),
  min_loci = 1,
  complete = FALSE,
  ...
) {
  rlang::check_dots_empty()
  # Check if the input is a gen_tibble
  stopifnot_gen_tibble(.x)
  type <- match.arg(type)

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

  # sum the quantities needed to compute tajimas d in each window, for each
  # group, in a single pass over the genotypes
  term_sums <- .tajimas_d_sums_by_group(.x,
    locus_first = win$locus_first,
    locus_last = win$locus_last,
    n_windows = nrow(win$windows)
  )

  res <- list()
  # now we can loop around each population to compute the windows
  for (i_grp in seq_len(ncol(term_sums$k_hat))) {
    grp_sums <- .tajimas_d_group_sums(term_sums, i_grp)
    # the number of loci with at least one valid allele copy in each window;
    # windows without such loci have no sums (as in windows_stats_generic())
    n_loci <- grp_sums[, "n_loci"]
    grp_sums[n_loci == 0, ] <- NA
    # compute tajimas d for each window from the sums
    stat <- tajimas_d_from_term_sums(grp_sums)
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

  groups <- .group_info(.x)
  names(res) <- groups$labels

  if (type == "matrix") {
    # one row per window: as the number of loci with data can differ between
    # groups (e.g. with an empty group), which would split a window over
    # several rows, n_loci is the number of loci in the window (the number
    # with data in each group is in the "tidy" and "list" outputs)
    n_loci <- as.numeric(win$hi - win$lo)
    n_loci[empty] <- NA
    if (complete) {
      n_loci[!win$complete] <- NA
    }
    res <- lapply(res, function(grp_res) {
      grp_res$n_loci <- n_loci
      grp_res
    })
    res <- bind_rows(res, .id = "group")
    res <-
      res %>% tidyr::pivot_wider(names_from = "group", values_from = "stat")
    return(res)
  } else if (type == "tidy") {
    res <- .bind_group_rows(res, groups)
    return(res)
  } else if (type == "list") {
    return(res)
  }
}

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
windows_pop_tajimas_d.data.frame <- function(.x, ...) {
  windows_pop_tajimas_d(restore_genotypes_tbl(.x, "windows_pop_tajimas_d"), ...)
}
