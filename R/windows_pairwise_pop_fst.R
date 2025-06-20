#' Compute pairwise Fst for a sliding window
#'
#' @description This function computes pairwise Fst for a sliding window across
#'   each chromosome.
#'
#' @param .x a grouped `gen_tibble` object
#' @param type type of object to return. One of "matrix" or "tidy". Default is
#'   "matrix". "matrix" returns a dataframe where each row is a window, followed
#'   by columns of Fst values for each pairwise population a and b comparison.
#'   "tidy" returns a tidy tibble of the same data in 'long' format, where each
#'   row is one window for one pairwise population a and b comparison.
#' @param method the method to use for calculating Fst, one of 'Hudson',
#'   'Nei87', and 'WC84'. See [pairwise_pop_fst()] for details.
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
#' - `fst_a.b`: the pairwise Fst value for the population a and b (there will be
#'   multiple such columns if there are more than two populations) or a tidy
#'   tibble with the following columns:
#' - `chromosome`: the chromosome for the window
#' - `start`: the starting locus of the window
#' - `end`: the ending locus of the window
#' - `stat_name`: the name of population a and b used in the pairwise Fst
#'   calculation (e.g. "fst_pop1.pop2")
#' - `value`: the pairwise Fst value for the population a and b
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>%
#'   group_by(population) %>%
#'   windows_pairwise_pop_fst(
#'     window_size = 3, step_size = 2,
#'     size_unit = "snp", min_loci = 2
#'   )
#'
windows_pairwise_pop_fst <- function(.x,
                                     type = c("matrix", "tidy"),
                                     method = c("Hudson", "Nei87", "WC84"),
                                     window_size,
                                     step_size,
                                     size_unit = c("snp", "bp"),
                                     min_loci = 1,
                                     complete = FALSE) {
  # Check if the input is a gen_tibble
  stopifnot_gen_tibble(.x)
  method <- match.arg(method)
  type <- match.arg(type)

  # create the pairwise Fst by locus, saving numerator and denominator
  pair_fst <- pairwise_pop_fst(.x,
    method = "Hudson", return_num_dem = TRUE,
    by_locus = TRUE
  )

  # for each column (i.e. each statistics), run the window analysis
  compn <- colnames(pair_fst$Fst_by_locus_num)
  res <- NULL
  for (i_col in seq_along(compn)) {
    # get the column name
    col_name <- compn[i_col]
    # compute the windows for numerator
    window_num <- windows_stats_generic(
      .x = pair_fst$Fst_by_locus_num[, i_col],
      loci_table = show_loci(.x),
      window_size = window_size,
      step_size = step_size,
      size_unit = size_unit,
      min_loci = min_loci,
      complete = complete
    )
    # same for the denominator
    window_dem <- windows_stats_generic(
      .x = pair_fst$Fst_by_locus_den[, i_col],
      loci_table = show_loci(.x),
      window_size = window_size,
      step_size = step_size,
      size_unit = size_unit,
      min_loci = min_loci,
      complete = complete
    )
    # compute the Fst for the window
    window_fst <- data.frame(stat = window_num$stat / window_dem$stat)
    names(window_fst) <- paste0("fst_", col_name)
    # if res is null
    if (is.null(res)) {
      # create the data frame
      res <- dplyr::bind_cols(
        window_num[, 1:3],
        window_fst
      )
    } else {
      res <- dplyr::bind_cols(res, window_fst)
    }
  }
  if (type == "matrix") {
    return(res)
  } else if (type == "tidy") {
    cols <- names(res)[names(res) != c("chromosome", "start", "end")]
    res <-
      res %>% tidyr::pivot_longer(cols = all_of(cols), names_to = "stat_name")
    return(res)
  }
}
