#' Create a Quality Control report for individuals
#'
#'#' Return QC information to assess loci (Observed heterozygosity and missingness).
#'
#' @param .x a [`gen_tibble`] object.
#' @param cutoff an optional numeric, a threshold of relatedness for the sample
#' @param ... further arguments to pass
#' @returns a tibble with 2 elements: het_obs and missingness
#' @export


indiv_qc_report <- function(.x, cutoff = NULL, ...){
  rlang::check_dots_empty()

  if(is.null(cutoff)){
    indiv_qc_report <- .x %>% reframe(het_obs = indiv_het_obs(.x),
                                      missingness = indiv_missingness(.x,as_counts=FALSE))
  } else {

    king <- gt_king(.x, as_matrix = TRUE)

    relatives <- filter_high_relatedness(matrix = king, .x = .x, cutoff = cutoff,...)

    indiv_qc_report <- .x %>% reframe(het_obs = indiv_het_obs(.x),
                                      missingness = indiv_missingness(.x,as_counts=FALSE))
    indiv_qc_report$to_keep <- relatives[[3]]
    indiv_qc_report$id <- .x$id
  }

  class(indiv_qc_report) <- c("indiv_qc_report",class(indiv_qc_report))
  indiv_qc_report
}


#' @export
autoplot.indiv_qc_report <- function(object, type = c("scatter"),miss_threshold = NULL, ...){

  rlang::check_dots_empty()

  miss_threshold <- if(is.null(miss_threshold)){
    0.05
  } else {
    miss_threshold
  }


  type <- match.arg(type)

  if (type == "scatter") {
    final_plot <- autoplot_indiv_qc_report(object,miss_threshold)
  } else {
    stop("Invalid type argument. Please choose from 'scatter' or 'relatedness'")
  }

  return(final_plot)

}

autoplot_indiv_qc_report <- function(object, miss_threshold = miss_threshold){

  miss_threshold <- miss_threshold

  mean_val <- mean(object$het_obs)
  sd_val <- stats::sd(object$het_obs)

  upper <- mean_val + 3*(sd_val)
  lower <- mean_val - 3*(sd_val)

  mid_upper <- mean_val + 2*(sd_val)
  mid_lower <- mean_val - 2*(sd_val)

  final_plot <- ggplot2::ggplot(object,ggplot2::aes(x=.data$missingness,y=.data$het_obs))+ggplot2::geom_point()+ggplot2::labs(x="Missingness",y="Observed Heterozygosity",title="Heterozygosity and missingness by individual")+ ggplot2::geom_vline(xintercept= miss_threshold, lty=2, col="red")+ggplot2::geom_hline(yintercept = c(upper,lower,mid_upper,mid_lower),lty=2,col="blue")

}
