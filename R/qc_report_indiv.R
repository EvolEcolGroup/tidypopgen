#' Create a Quality Control report for individuals
#'
#'#' Return QC information to assess loci (Observed heterozygosity and missingness).
#'
#' @param .x a [`gen_tibble`] object.
#' @param kings_threshold an optional numeric, a threshold of relatedness for the sample
#' @param ... further arguments to pass
#' @returns a tibble with 2 elements: het_obs and missingness
#' @export


qc_report_indiv <- function(.x, kings_threshold = NULL, ...){
  rlang::check_dots_empty()

  if(is.null(kings_threshold)){
    qc_report_indiv <- .x %>% reframe(het_obs = indiv_het_obs(.x),
                                      missingness = indiv_missingness(.x,as_counts=FALSE))
  } else {

    king <- gt_king(.x, as_matrix = TRUE)

    relatives <- filter_high_relatedness(matrix = king, .x = .x, kings_threshold = kings_threshold,...)

    qc_report_indiv <- .x %>% reframe(het_obs = indiv_het_obs(.x),
                                      missingness = indiv_missingness(.x,as_counts=FALSE))
    qc_report_indiv$to_keep <- relatives[[3]]
    qc_report_indiv$id <- .x$id
    attr(qc_report_indiv$to_keep, "king") <- king

  }

  class(qc_report_indiv) <- c("qc_report_indiv",class(qc_report_indiv))
  qc_report_indiv
}


#' Autoplots for `qc_report_indiv` objects
#'
#' For `qc_report_indiv`, the following types of plots are available:
#' - `scatter`: a plot of missingness and observed heterozygosity within
#' individuals.
#' - `relatedness`: a histogram of paired kinship coefficients
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are
#' not customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' @param object an object of class `qc_report_indiv`
#' @param type the type of plot (`scatter`,`relatedness`)
#' @param miss_threshold a threshold for the accepted rate of missingness within
#' individuals
#' @param kings_threshold an optional numeric, a threshold of relatedness for the sample
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @export
autoplot.qc_report_indiv <- function(object, type = c("scatter", "relatedness"),miss_threshold = NULL, kings_threshold = kings_threshold, ...){

  rlang::check_dots_empty()

  miss_threshold <- if(is.null(miss_threshold)){
    0.05
  } else {
    miss_threshold
  }

  type <- match.arg(type)

  if (type == "scatter") {
    final_plot <- autoplot_qc_report_indiv(object,miss_threshold)
  } else if (type == "relatedness") {
    final_plot <- autoplot_qc_report_indiv_king(object,kings_threshold)
  } else {
    stop("Invalid type argument. Please choose from 'scatter' or 'relatedness'")
  }

  return(final_plot)

}

autoplot_qc_report_indiv <- function(object, miss_threshold = miss_threshold){

  miss_threshold <- miss_threshold

  mean_val <- mean(object$het_obs)
  sd_val <- stats::sd(object$het_obs)

  upper <- mean_val + 3*(sd_val)
  lower <- mean_val - 3*(sd_val)

  mid_upper <- mean_val + 2*(sd_val)
  mid_lower <- mean_val - 2*(sd_val)

  final_plot <- ggplot2::ggplot(object,ggplot2::aes(x=.data$missingness,y=.data$het_obs))+ggplot2::geom_point()+ggplot2::labs(x="Missingness",y="Observed Heterozygosity",title="Heterozygosity and missingness by individual")+ ggplot2::geom_vline(xintercept= miss_threshold, lty=2, col="red")+ggplot2::geom_hline(yintercept = c(upper,lower,mid_upper,mid_lower),lty=2,col="blue")

}

autoplot_qc_report_indiv_king <- function(object, kings_threshold = kings_threshold){

  king <- as.data.frame(attr(object$to_keep, "king"))
  num_samples <- nrow(king)
  king$row <- colnames(king)

  #format into 3 columns: ID1, ID2, and their relatedness coefficient
  king <- tidyr::pivot_longer(king, cols = !row, names_to = "Column",
                              values_to = "Value")
  king <- dplyr::filter(king, king$row < king$Column)

  colnames(king) <- c("ID1","ID2","kinship")

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
  total_pairs <- num_samples * (num_samples -1)/2

  if (total_pairs != nrow(king_sorted)){
    stop("Relatedness matrix must be symmetric ")
  }

  p <- ggplot2::ggplot(king_sorted, ggplot2::aes(x=.data$kinship)) +
    ggplot2::geom_histogram(bins = 40) +
    ggplot2::labs(x="KING robust kinship estimator",y = "Number of pairs",
                  title = "Distribution of paired kinship coefficients") +
    ggplot2::geom_vline(xintercept = kings_threshold, lty=2, col="red")


  }

