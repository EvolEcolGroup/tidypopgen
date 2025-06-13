#' Create a Quality Control report for individuals
#'
#' Return QC information to assess loci (Observed heterozygosity and
#' missingness).
#'
#' Providing the parameter kings_threshold will return two additional columns,
#' 'id' containing the ID of individuals, and 'to_keep' a logical vector
#' describing whether the individual should be removed to retain the largest
#' possible set of individuals with no relationships above the threshold. The
#' calculated pairwise KING relationship matrix is also returned as an attribute
#' of 'to_keep'. The kings_threshold parameter can be either a numeric KING
#' kinship coefficient or a string of either "first" or "second", to remove any
#' first degree or second degree relationships from the dataset. This second
#' option is similar to using  --unrelated --degree 1 or --unrelated --degree 2
#' in KING.
#'
#' @param .x either a [`gen_tibble`] object or a grouped [`gen_tibble`] (as
#'   obtained by using [dplyr::group_by()])
#' @param kings_threshold an optional numeric giving a KING kinship coefficient,
#' or one of:
#'   - "first": removing first degree relatives, equivalent to a kinship
#'   coefficient of 0.177 or more
#'   - "second": removing second degree relatives, equivalent to a kinship
#'   coefficient of 0.088 or more
#' @param ... further arguments to pass
#' @returns If no kings_threshold is provided, a tibble with 2 elements: het_obs
#'   and missingness. If kings_threshold is provided, a tibble with 4 elements:
#'   het_obs, missingness, id and to_keep.
#' @rdname qc_report_indiv
#' @export
#' @examples
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' example_gt <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Get QC report for individuals
#' example_gt %>% qc_report_indiv()
#'
#' # Get QC report with kinship filtering
#' example_gt %>% qc_report_indiv(kings_threshold = "first")
qc_report_indiv <- function(.x, ...) {
  UseMethod("qc_report_indiv", .x)
}


#' @export
#' @rdname qc_report_indiv
qc_report_indiv.tbl_df <- function(.x, kings_threshold = NULL, ...) {
  rlang::check_dots_empty()

  if (!is.null(kings_threshold)) {
    if (kings_threshold %in% c("first", "second")) {
      if (kings_threshold == "first") {
        kings_threshold <- 0.177
      } else if (kings_threshold == "second") {
        kings_threshold <- 0.088
      }
    } else if (!is.numeric(kings_threshold)) {
      stop("kings_threshold must be a numeric or one of 'first' or 'second'")
    }
  }

  n_loci <- nrow(show_loci(.x))
  qc_report <- .x %>%
    indiv_het_obs(, as_counts = TRUE) %>%
    as_tibble() %>%
    reframe(
      het_obs = .data$het_n / (n_loci - .data$na_n),
      missingness = .data$na_n / n_loci
    )
  if (!is.null(kings_threshold)) {
    # calculate the kinship matrix
    king <- pairwise_king(.x, as_matrix = TRUE)

    relatives <- filter_high_relatedness(
      matrix = king,
      .x = .x,
      kings_threshold = kings_threshold,
      ...
    )

    qc_report$to_keep <- relatives[[3]]
    qc_report$id <- .x$id
    attr(qc_report$to_keep, "king") <- king
  }
  class(qc_report) <- c("qc_report_indiv", class(qc_report))
  qc_report
}


#' @export
#' @rdname qc_report_indiv
qc_report_indiv.grouped_df <- function(.x, kings_threshold = NULL, ...) {
  rlang::check_dots_empty()

  if (!is.null(kings_threshold)) {
    if (is.numeric(kings_threshold)) {
      kings_threshold <- kings_threshold
    } else if (!is.numeric(kings_threshold) && !kings_threshold %in% c("first", "second")) { # nolint
      stop("kings_threshold must be a numeric or one of 'first' or 'second'")
    } else {
      kings_threshold <- match.arg(kings_threshold, c("first", "second"))
    }

    if (kings_threshold == "first") {
      kings_threshold <- 0.177
    } else if (kings_threshold == "second") {
      kings_threshold <- 0.088
    }
  }

  n_loci <- nrow(show_loci(.x))
  qc_report <- .x %>%
    # TODO do we need this???
    ungroup() %>%
    indiv_het_obs(, as_counts = TRUE) %>%
    as_tibble() %>%
    reframe(
      het_obs = .data$het_n / (n_loci - .data$na_n),
      missingness = .data$na_n / n_loci
    ) %>%
    mutate(id = .x$id)
  # re-group according to original .x grouping
  qc_report$group <- .x[[dplyr::group_vars(.x)]]

  if (!is.null(kings_threshold)) {
    # find grouping levels
    grouping_order <- group_keys(.x)

    # calculate the kinship matrix for each population
    king <- .x %>% group_map(~ pairwise_king(.x, as_matrix = TRUE))
    names(king) <- as.matrix(grouping_order)[, 1]

    # create relatives, an empty list
    relatives <- vector("list")

    # remove individuals using filter_high_relatedness for each population
    indivs_to_keep <- function(relatives, king) {
      for (i in seq_along(king)) {
        output <- filter_high_relatedness(
          king[[i]],
          kings_threshold = kings_threshold
        )
        # add the output to the list
        relatives[[i]] <- output[[3]]
        names(relatives)[i] <- names(king)[i]
      }
      return(relatives) # nolint
    }

    relatives <- indivs_to_keep(relatives, king)

    # reorder relatives
    relatives <- relatives[unique(qc_report$group)]
    all_flags <- unlist(relatives)
    # add to_remove
    qc_report$to_keep <- as.vector(all_flags)
    attr(qc_report$to_keep, "king") <- king
  }

  class(qc_report) <- c("qc_report_indiv", class(qc_report))
  return(qc_report)
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
#' @examples
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' example_gt <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Create QC report for individuals
#' indiv_report <- example_gt %>% qc_report_indiv()
#'
#' # Autoplot missingness and observed heterozygosity
#' autoplot(indiv_report, type = "scatter", miss_threshold = 0.1)
#'
#' # Create QC report with kinship filtering
#' indiv_report_rel <-
#'   example_gt %>% qc_report_indiv(kings_threshold = "second")
#'
#' # Autoplot relatedness
#' autoplot(indiv_report_rel, type = "relatedness", kings_threshold = "second")
#'
autoplot.qc_report_indiv <- function(
    object,
    type = c("scatter", "relatedness"),
    miss_threshold = 0.05,
    kings_threshold = NULL,
    ...) {
  rlang::check_dots_empty()

  if (!is.null(kings_threshold)) {
    if (is.numeric(kings_threshold)) {
      kings_threshold <- kings_threshold
    } else if (!is.numeric(kings_threshold) && !kings_threshold %in% c("first", "second")) { # nolint
      stop("kings_threshold must be a numeric or one of 'first' or 'second'")
    } else {
      kings_threshold <- match.arg(kings_threshold, c("first", "second"))
    }


    if (kings_threshold == "first") {
      kings_threshold <- 0.177
    } else if (kings_threshold == "second") {
      kings_threshold <- 0.088
    }
  }

  type <- match.arg(type)

  if (type == "scatter") {
    report_plot <- autoplot_qc_report_indiv(
      object,
      miss_threshold = miss_threshold
    )
  } else if (type == "relatedness") {
    report_plot <- autoplot_qc_report_indiv_king(
      object,
      kings_threshold = kings_threshold
    )
  }
  return(report_plot)
}

autoplot_qc_report_indiv <- function(object, miss_threshold) {
  mean_val <- mean(object$het_obs)
  sd_val <- stats::sd(object$het_obs)

  upper <- mean_val + 3 * (sd_val)
  lower <- mean_val - 3 * (sd_val)

  mid_upper <- mean_val + 2 * (sd_val)
  mid_lower <- mean_val - 2 * (sd_val)

  report_plot <- ggplot2::ggplot(
    object,
    ggplot2::aes(
      x = .data$missingness,
      y = .data$het_obs
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Missingness",
      y = "Observed Heterozygosity",
      title = "Heterozygosity and missingness by individual"
    ) +
    ggplot2::geom_vline(xintercept = miss_threshold, lty = 2, col = "red") +
    ggplot2::geom_hline(
      yintercept = c(upper, lower, mid_upper, mid_lower),
      lty = 2,
      col = "blue"
    )
  return(report_plot)
}

autoplot_qc_report_indiv_king <- function(
    object,
    kings_threshold) {
  if (
    inherits(attr(object$to_keep, "king"), "matrix") ||
      inherits(attr(object$to_keep, "king"), "array")
  ) {
    p <- king_hist_plot(attr(object$to_keep, "king"), kings_threshold)
  } else if (inherits(attr(object$to_keep, "king"), "list")) {
    # by group

    # apply plotting function to each kings matrix in the list
    p <- lapply(
      FUN = king_hist_plot, X = attr(object$to_keep, "king"),
      kings_threshold = kings_threshold
    )
  }

  return(p)
}

king_hist_plot <- function(object, kings_threshold) {
  king <- as.data.frame(object)
  num_samples <- nrow(king)
  king$row <- colnames(king)

  # format into 3 columns: ID1, ID2, and their relatedness coefficient
  king <- tidyr::pivot_longer(
    king,
    cols = !row,
    names_to = "Column",
    values_to = "Value"
  )
  king <- dplyr::filter(king, king$row < king$Column)
  colnames(king) <- c("ID1", "ID2", "kinship")

  # remove duplication from the new df
  king_sorted <- king %>%
    mutate(
      row_min = pmin(.data$ID1, .data$ID2),
      row_max = pmax(.data$ID1, .data$ID2)
    ) %>%
    dplyr::select(-dplyr::all_of(c("ID1", "ID2"))) %>%
    distinct(.data$row_min, .data$row_max, .keep_all = TRUE) %>%
    rename("ID1" = "row_min", "ID2" = "row_max")

  # remove cases of individuals relatedness with themselves
  king_sorted <- king_sorted %>% filter(.data$ID1 != .data$ID2)

  # add a check for correct number of pairs
  total_pairs <- num_samples * (num_samples - 1) / 2

  if (total_pairs != nrow(king_sorted)) {
    stop("Relatedness matrix must be symmetric ")
  }

  p <- ggplot2::ggplot(king_sorted, ggplot2::aes(x = .data$kinship)) +
    ggplot2::geom_histogram(bins = 40) +
    ggplot2::labs(
      x = "KING robust kinship estimator",
      y = "Number of pairs",
      title = "Distribution of paired kinship coefficients"
    ) +
    ggplot2::geom_vline(xintercept = kings_threshold, lty = 2, col = "red")
  return(p)
}
