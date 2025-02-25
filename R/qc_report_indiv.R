#'Create a Quality Control report for individuals
#'
#'Return QC information to assess loci (Observed heterozygosity and
#'missingness).
#'
#'@param .x either a [`gen_tibble`] object or a grouped [`gen_tibble`] (as
#'  obtained by using [dplyr::group_by()])
#'@param kings_threshold an optional numeric, a threshold of relatedness for the
#'  sample
#'@param ... further arguments to pass
#'@returns a tibble with 2 elements: het_obs and missingness
#'@rdname qc_report_indiv
#'@export
qc_report_indiv <- function(.x, ...) {
  UseMethod("qc_report_indiv", .x)
}


#' @export
#' @rdname qc_report_indiv
qc_report_indiv.tbl_df <- function(.x, kings_threshold = NULL, ...) {
  rlang::check_dots_empty()

  if (is.null(kings_threshold)) {
    qc_report_indiv <- .x %>%
      reframe(het_obs = indiv_het_obs(.x),
              missingness = indiv_missingness(.x, as_counts = FALSE))
  } else {

    # calculate the kinship matrix
    king <- pairwise_king(.x, as_matrix = TRUE)

    relatives <- filter_high_relatedness(matrix = king, .x = .x,
                                         kings_threshold = kings_threshold, ...)

    qc_report_indiv <- .x %>%
      reframe(het_obs = indiv_het_obs(.x),
              missingness = indiv_missingness(.x, as_counts = FALSE))
    qc_report_indiv$to_keep <- relatives[[3]]
    qc_report_indiv$id <- .x$id
    attr(qc_report_indiv$to_keep, "king") <- king

  }
  class(qc_report_indiv) <- c("qc_report_indiv", class(qc_report_indiv))
  qc_report_indiv
}


#' @export
#' @rdname qc_report_indiv
qc_report_indiv.grouped_df <- function(.x, kings_threshold = NULL, ...) {
  rlang::check_dots_empty()

  if (is.null(kings_threshold)) {
    qc_report_indiv <- .x %>%
      ungroup() %>%
      reframe(het_obs = indiv_het_obs(.x),
              missingness = indiv_missingness(.x, as_counts = FALSE))
  } else {

  # find grouping levels
  grouping_order <- group_keys(.x)

  # calculate the kinship matrix for each population
  king <- .x %>% group_map(~pairwise_king(.x, as_matrix = TRUE))
  names(king) <- as.matrix(grouping_order)[, 1]

  # create relatives, an empty list
  relatives <- vector("list")

  # remove individuals using filter_high_relatedness for each population
  indivs_to_keep <- function(relatives, king) {
    for (i in 1:length(king)){
      output <- filter_high_relatedness(king[[i]],
                                        kings_threshold = kings_threshold)
      # add the output to the list
      relatives[[i]] <- output[[3]]
      names(relatives)[i] <- names(king)[i]
    }
    return(relatives)
  }

  relatives <- indivs_to_keep(relatives, king)

  # create qc_report_indiv
  qc_report_indiv <- .x %>%
    ungroup() %>%
    reframe(het_obs = indiv_het_obs(.x),
            missingness = indiv_missingness(.x, as_counts = FALSE)) %>%
    mutate(id = .x$id)
  # re-group according to original .x grouping
  qc_report_indiv$group <- .x[[dplyr::group_vars(.x)]]
  # reorder relatives
  relatives <- relatives[unique(qc_report_indiv$group)]
  all_flags <- unlist(relatives)
  # add to_remove
  qc_report_indiv$to_keep <- as.vector(all_flags)
  attr(qc_report_indiv$to_keep, "king") <- king

  }

  class(qc_report_indiv) <- c("qc_report_indiv", class(qc_report_indiv))
  return(qc_report_indiv)

}




#' Autoplots for `qc_report_indiv` objects
#'
#' For `qc_report_indiv`, the following types of plots are available:
#' - `scatter`: a plot of missingness and observed heterozygosity within
#' individuals.
#' - `relatedness`: a histogram of paired kinship coefficients
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are not
#' customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' @param object an object of class `qc_report_indiv`
#' @param type the type of plot (`scatter`,`relatedness`)
#' @param miss_threshold a threshold for the accepted rate of missingness within
#'   individuals
#' @param kings_threshold an optional numeric, a threshold of relatedness for
#'   the sample
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @export
autoplot.qc_report_indiv <- function(object, type = c("scatter", "relatedness"),
                                     miss_threshold = NULL,
                                     kings_threshold = kings_threshold, ...) {

  rlang::check_dots_empty()

  miss_threshold <- if (is.null(miss_threshold)) {
    0.05
  } else {
    miss_threshold
  }

  type <- match.arg(type)

  if (type == "scatter") {
    final_plot <- autoplot_qc_report_indiv(object,
                                           miss_threshold = miss_threshold)
  } else if (type == "relatedness") {
    final_plot <- autoplot_qc_report_indiv_king(object,
                                                kings_threshold = kings_threshold)
  }
  return(final_plot)
}

autoplot_qc_report_indiv <- function(object, miss_threshold = miss_threshold) {

  miss_threshold <- miss_threshold

  mean_val <- mean(object$het_obs)
  sd_val <- stats::sd(object$het_obs)

  upper <- mean_val + 3 * (sd_val)
  lower <- mean_val - 3 * (sd_val)

  mid_upper <- mean_val + 2 * (sd_val)
  mid_lower <- mean_val - 2 * (sd_val)

  final_plot <- ggplot2::ggplot(object, ggplot2::aes(x = .data$missingness,
                                                     y = .data$het_obs)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Missingness", y = "Observed Heterozygosity",
                  title = "Heterozygosity and missingness by individual") +
    ggplot2::geom_vline(xintercept = miss_threshold, lty = 2, col = "red") +
    ggplot2::geom_hline(yintercept = c(upper, lower, mid_upper, mid_lower),
                        lty = 2, col = "blue")
  return(final_plot)
}

autoplot_qc_report_indiv_king <- function(object, kings_threshold = kings_threshold) {

  if (inherits(attr(object$to_keep, "king"), "matrix") || inherits(attr(object$to_keep, "king"), "array")) {

    king <- as.data.frame(attr(object$to_keep, "king"))
    num_samples <- nrow(king)
    king$row <- colnames(king)

    #format into 3 columns: ID1, ID2, and their relatedness coefficient
    king <- tidyr::pivot_longer(king, cols = !row, names_to = "Column",
                                values_to = "Value")
    king <- dplyr::filter(king, king$row < king$Column)
    colnames(king) <- c("ID1", "ID2", "kinship")

    #remove duplication from the new df
    king_sorted <- king %>%
      mutate(
        row_min = pmin(.data$ID1, .data$ID2),
        row_max = pmax(.data$ID1, .data$ID2)
      ) %>%
      dplyr::select(-dplyr::all_of(c("ID1", "ID2"))) %>%
      distinct(.data$row_min, .data$row_max, .keep_all = TRUE) %>%
      rename("ID1" = .data$row_min, "ID2" = .data$row_max)

    # remove cases of individuals relatedness with themselves
    king_sorted <- king_sorted %>% filter(.data$ID1 != .data$ID2)

    #add a check for correct number of pairs
    total_pairs <- num_samples * (num_samples - 1) / 2

    if (total_pairs != nrow(king_sorted)) {
      stop("Relatedness matrix must be symmetric ")
    }

    p <- ggplot2::ggplot(king_sorted, ggplot2::aes(x = .data$kinship)) +
      ggplot2::geom_histogram(bins = 40) +
      ggplot2::labs(x = "KING robust kinship estimator", y = "Number of pairs",
                    title = "Distribution of paired kinship coefficients") +
      ggplot2::geom_vline(xintercept = kings_threshold, lty = 2, col = "red")


  } else if (inherits(attr(object$to_keep, "king"), "list")) { # by group

    # create a plotting function as above
    hist_plot <- function(object) {
      king <- as.data.frame(object)
      num_samples <- nrow(king)
      king$row <- colnames(king)
      king <- tidyr::pivot_longer(king, cols = !row, names_to = "Column",
                                  values_to = "Value")
      king <- dplyr::filter(king, king$row < king$Column)
      colnames(king) <- c("ID1", "ID2", "kinship")
      king_sorted <- king %>%
        mutate(
          row_min = pmin(.data$ID1, .data$ID2),
          row_max = pmax(.data$ID1, .data$ID2)
        ) %>%
        dplyr::select(-dplyr::all_of(c("ID1", "ID2"))) %>%
        distinct(.data$row_min, .data$row_max, .keep_all = TRUE) %>%
        rename("ID1" = .data$row_min, "ID2" = .data$row_max)
      king_sorted <- king_sorted %>% filter(.data$ID1 != .data$ID2)
      total_pairs <- num_samples * (num_samples - 1) / 2
      if (total_pairs != nrow(king_sorted)) {
        stop("Relatedness matrix must be symmetric ")
      }
      p <- ggplot2::ggplot(king_sorted, ggplot2::aes(x = .data$kinship)) +
        ggplot2::geom_histogram(bins = 40) +
        ggplot2::labs(x = "KING robust kinship estimator",
                      y = "Number of pairs",
                      title = "Distribution of paired kinship coefficients") +
        ggplot2::geom_vline(xintercept = kings_threshold, lty = 2, col = "red")

      return(p)
    }

    # apply plotting function to each kings matrix in the list
    p <- lapply(FUN = hist_plot, X = attr(object$to_keep, "king"))

  }

  return(p)
}
