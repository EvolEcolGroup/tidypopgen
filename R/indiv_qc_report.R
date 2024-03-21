#' Create a Quality Control report for individuals
#'
#'#' Return QC information to assess loci (Observed heterozygosity and missingness).
#'
#' @param .x a [`gen_tibble`] object.
#' @returns a tibble with 2 elements: het_obs and missingness
#' @export
indiv_qc_report <- function(.x,...){
  indiv_qc_report <- .x %>% reframe(het_obs = indiv_het_obs(.x),
                              missingness = indiv_missingness(.x,as_counts=FALSE))
  class(indiv_qc_report) <- c("indiv_qc_report",class(indiv_qc_report))
  indiv_qc_report
}

#' @export
autoplot.indiv_qc_report <- function(.x, type = c("scatter"),miss_threshold = NULL,...){

  miss_threshold <- if(is.null(miss_threshold)){
    0.05
  } else {
    miss_threshold
  }

  type <- match.arg(type)

  if (type == "scatter") {
    final_plot <- autoplot_indiv_qc_report(.x,miss_threshold)
  } else {
    stop("Invalid type argument. Please choose 'scatter'")
  }

  return(final_plot)

}


autoplot_indiv_qc_report <- function(.x, miss_threshold = miss_threshold,...){

  miss_threshold <- miss_threshold

  mean_val <- mean(.x$het_obs)
  sd_val <- sd(.x$het_obs)

  upper <- mean_val + 3*(sd_val)
  lower <- mean_val - 3*(sd_val)

  mid_upper <- mean_val + 2*(sd_val)
  mid_lower <- mean_val - 2*(sd_val)

  final_plot <- ggplot(.x,aes(x=missingness,y=het_obs))+geom_point()+labs(x="Missingness",y="Observed Heterozygosity",title="Heterozygosity and missingness by individual")+ geom_vline(xintercept= miss_threshold, lty=2, col="red")+geom_hline(yintercept = c(upper,lower,mid_upper,mid_lower),lty=2,col="blue")

}
