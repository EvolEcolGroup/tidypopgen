#' Create a Quality Control report for loci
#'
#' Return QC information to assess loci (MAF, missingness and HWE test).
#'
#' @param .x a [`gen_tibble`] object.
#' @param ... further arguments to pass to [HardyWeinberg::HWExact()] for changing
#' the HWE test.
#' @returns a tibble with 3 elements: maf, missingness and hwe_p
#' @export
loci_qc_report <- function (.x, ...){
  qc_report <- .x %>% reframe(snp_id = show_loci_names(.x),maf=loci_maf(.data$genotypes),
                              missingness = loci_missingness(.data$genotypes),
                    hwe_p = loci_hwe(.data$genotypes, ...))
 class(qc_report) <- c("loci_qc_report",class(qc_report))
 qc_report
}


#' @export
autoplot.loci_qc_report <- function(object, type = c("overview","all","missing","missing low maf","missing high maf","maf","hwe","significant hwe"), maf_threshold = NULL, miss_threshold = NULL, p_val = NULL,...) {

  type <- match.arg(type)

  rlang::check_dots_empty()

  maf_threshold <- if(is.null(maf_threshold)){
    0.05
  } else {
    maf_threshold
  }

  miss_threshold <- if(is.null(miss_threshold)){
    0.01
  } else {
    miss_threshold
  }

  p_val <- if(is.null(p_val)){
    0.01
  } else {
    p_val
  }

  logp <- -log10(p_val)

  if (type == "overview") {
    final_plot <- autoplot_l_qc_overview(object, maf_threshold, miss_threshold)
  } else if (type == "all") {
    final_plot <- autoplot_l_qc_all(object, maf_threshold, miss_threshold, p_val,logp)
  } else if (type == "missing") {
    final_plot <- autoplot_l_qc_missing(object, miss_threshold)
  } else if (type == "missing low maf") {
    final_plot <- autoplot_l_qc_missing_low_maf(object, maf_threshold, miss_threshold)
  } else if (type == "missing high maf") {
    final_plot <- autoplot_l_qc_missing_high_maf(object, maf_threshold, miss_threshold)
  } else if (type == "maf") {
    final_plot <- autoplot_l_qc_maf(object, maf_threshold)
  } else if (type == "hwe") {
    final_plot <- autoplot_l_qc_hwe(object,logp)
  } else if (type == "significant hwe") {
    final_plot <- autoplot_l_qc_sig_hwe(object,p_val,logp)
  } else {
    stop("Invalid type argument. Please choose from 'overview','all','maf','hwe','significant hwe'")
  }

  return(final_plot)
}





autoplot_l_qc_all <- function(object, maf_threshold = maf_threshold, miss_threshold = miss_threshold, p_val = p_val, logp = logp,...){

  qc_report <- object

  #Missingness (according to MAF thresholds)
  qc_highmaf <- subset(qc_report, qc_report$maf > maf_threshold)
  qc_lowmaf <- subset(qc_report, qc_report$maf <= maf_threshold)

  p_highMAF <- ggplot2::ggplot(qc_highmaf,ggplot2::aes(x=.data$missingness))+ggplot2::geom_histogram(position = "dodge",binwidth = 0.005,fill="#66C2A5") + ggplot2::labs(x="Proportion of missing data",y="Number of SNPs",title = paste("SNPs with MAF > ", maf_threshold))+
    ggplot2::geom_vline(xintercept=miss_threshold, lty=2, col="red")+
    ggplot2::scale_color_brewer(palette = "Dark2")

  p_lowMAF <- ggplot2::ggplot(qc_lowmaf,ggplot2::aes(x=.data$missingness))+ggplot2::geom_histogram(position = "dodge",binwidth=0.005,fill="#66C2A5") + ggplot2::labs(x="Proportion of missing data",y="Number of SNPs", title = paste("SNPs with MAF < ", maf_threshold)) +
    ggplot2::geom_vline(xintercept=miss_threshold, lty=2, col="red")+
    ggplot2::scale_color_brewer(palette = "Dark2")

  mafmiss <- patchwork::wrap_plots(p_highMAF + p_lowMAF)


  #Minor allele frequency distribution
  #maf <- ggplot2::ggplot(qc_report,ggplot2::aes(x=.data$maf))+ggplot2::geom_histogram(binwidth= 0.01,fill="#66C2A5")+ ggplot2::labs(x="Minor allele frequency",y="Number of SNPs", title = "Minor allele frequency distribution")+ ggplot2::geom_vline(xintercept = maf_threshold, lty=2, col="red")

  #Hardy weinberg exact test p-val distribution
  qc_report$hwe_p_log <- -log10(qc_report$hwe_p)
  qc_lowhwe <- subset(qc_report,qc_report$hwe_p < p_val)

  hwe_all <- ggplot2::ggplot(qc_report,ggplot2::aes(x=.data$hwe_p_log))+ggplot2::geom_histogram(binwidth = 0.5,fill="#66C2A5")+ ggplot2::labs(x=expression("-log"[10]*" of HWE exact p-value"),y="Number of SNPs", title = "Hardy-Weinberg exact test")+ ggplot2::geom_vline(xintercept= logp, lty=2, col="red")
  hwe_low <- ggplot2::ggplot(qc_lowhwe,ggplot2::aes(x=.data$hwe_p_log))+ggplot2::geom_histogram(binwidth = 0.5,fill="#66C2A5")+ ggplot2::labs(x=expression("-log"[10]* " of HWE exact p-value"),y="Number of SNPs", title = paste("HWE exact p-value <", p_val))+ ggplot2::geom_vline(xintercept= logp, lty=2, col="red")

  hwes <- patchwork::wrap_plots(hwe_all,hwe_low)

  final_plot_all <- mafmiss / hwes

}


autoplot_l_qc_overview <- function(object,maf_threshold=maf_threshold, miss_threshold = miss_threshold,...){

  qc_report <- object

  qc_lowmaf <- subset(qc_report, qc_report$maf <= maf_threshold)
  qc_lowhwe <- subset(qc_report,qc_report$hwe_p < 0.01)


  maf_fail <- c(qc_lowmaf$snp_id)
  hwe_fail <- c(qc_lowhwe$snp_id)

  qc_highmissing <- subset(qc_report,qc_report$missingness>=0.05)
  missing_fail <- c(qc_highmissing$snp_id)

  list <- list(MAF = maf_fail, HWE = hwe_fail, Missing = missing_fail)

  final_plot_overview <- UpSetR::upset(UpSetR::fromList(list),order.by = "freq",main.bar.color="#66C2A5", matrix.color="#66C2A5",
                                       sets.bar.color="#FC8D62")
}

autoplot_l_qc_maf <- function(object,maf_threshold=maf_threshold,...){

  qc_report <- object

  #Minor allele frequency distribution
  maf <- ggplot2::ggplot(qc_report,ggplot2::aes(x=.data$maf))+ggplot2::geom_histogram(binwidth= 0.01,fill="#66C2A5")+ ggplot2::labs(x="Minor allele frequency",y="Number of SNPs", title = "Minor allele frequency distribution")+ ggplot2::geom_vline(xintercept = maf_threshold, lty=2, col="red")

}

autoplot_l_qc_hwe <- function(object,logp,...){

  qc_report <- object

  #Hardy weinberg exact test p-val distribution
  qc_report$hwe_p_log <- -log10(qc_report$hwe_p)

  hwe_all <- ggplot2::ggplot(qc_report,ggplot2::aes(x=.data$hwe_p_log))+ggplot2::geom_histogram(binwidth = 0.5,fill="#66C2A5")+ ggplot2::labs(x=expression("-log"[10]*" of HWE exact p-value"),y="Number of SNPs", title = "Hardy-Weinberg exact test")+ ggplot2::geom_vline(xintercept= logp, lty=2, col="red")

}

autoplot_l_qc_sig_hwe <- function(object,p_val=p_val,logp,...){

  qc_report <- object

  #Hardy weinberg exact test p-val distribution
  qc_report$hwe_p_log <- -log10(qc_report$hwe_p)
  qc_lowhwe <- subset(qc_report,qc_report$hwe_p < p_val)

  hwe_low <- ggplot2::ggplot(qc_lowhwe,ggplot2::aes(x=.data$hwe_p_log))+ggplot2::geom_histogram(binwidth = 0.5,fill="#66C2A5")+ ggplot2::labs(x=expression("-log"[10]* " of HWE exact p-value"),y="Number of SNPs", title = paste("HWE exact p-value <", p_val))+ ggplot2::geom_vline(xintercept= logp, lty=2, col="red")
}

autoplot_l_qc_missing <- function(object, miss_threshold = miss_threshold, ...){

  qc_report <- object

  missing <- ggplot2::ggplot(qc_report,ggplot2::aes(x=.data$missingness))+ggplot2::geom_histogram(position = "dodge",binwidth=0.005,fill="#66C2A5") + ggplot2::labs(x="Proportion of missing data",y="Number of SNPs", title = "SNP missingness") +
    ggplot2::geom_vline(xintercept=miss_threshold, lty=2, col="red")+
    ggplot2::scale_color_brewer(palette = "Dark2")

}


autoplot_l_qc_missing_low_maf <- function(object,maf_threshold=maf_threshold,miss_threshold = miss_threshold,...){

  qc_report <- object

  qc_lowmaf <- subset(qc_report, qc_report$maf <= maf_threshold)

  p_highMAF <- ggplot2::ggplot(qc_lowmaf,ggplot2::aes(x=.data$missingness))+ggplot2::geom_histogram(position = "dodge",binwidth=0.005,fill="#66C2A5") + ggplot2::labs(x="Proportion of missing data",y="Number of SNPs", title = paste("SNPs with MAF <",maf_threshold)) +
    ggplot2::geom_vline(xintercept=miss_threshold, lty=2, col="red")+
    ggplot2::scale_color_brewer(palette = "Dark2")
}

autoplot_l_qc_missing_high_maf <- function(object,maf_threshold=maf_threshold,miss_threshold = miss_threshold,...){

  qc_report <- object

  qc_highmaf <- subset(qc_report, qc_report$maf > maf_threshold)

  p_highMAF <- ggplot2::ggplot(qc_highmaf,ggplot2::aes(x=.data$missingness))+ggplot2::geom_histogram(position = "dodge",binwidth = 0.005,fill="#66C2A5") + ggplot2::labs(x="Proportion of missing data",y="Number of SNPs",title = paste("SNPs with MAF >",maf_threshold))+
    ggplot2::geom_vline(xintercept=miss_threshold, lty=2, col="red")+
    ggplot2::scale_color_brewer(palette = "Dark2")
}
