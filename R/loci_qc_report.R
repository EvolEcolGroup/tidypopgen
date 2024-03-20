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
  qc_report <- .x %>% reframe(snp_id = show_loci_names(.x),maf=loci_freq(.data$genotypes),
                              missingness = loci_missingness(.data$genotypes),
                    hwe_p = loci_hwe(.data$genotypes, ...))
 class(qc_report) <- c("loci_qc_report",class(qc_report))
 qc_report
}


#' @export
autoplot.loci_qc_report <- function(.x, type = c("overview","all"), ...) {
  type <- match.arg(type)


  if (type == "overview") {
    final_plot <- autoplot_l_qc_overview(.x)
  } else if (type == "all") {
    final_plot <- autoplot_l_qc_all(.x)
  } else {
    stop("Invalid type argument")
  }

  return(final_plot)
}





autoplot_l_qc_all <- function(.x,...){

  qc_report <- .x

  #Missingness (according to MAF thresholds)
  qc_highmaf <- subset(qc_report, qc_report$maf > 0.05)
  qc_lowmaf <- subset(qc_report, qc_report$maf <= 0.05)

  p_highMAF <- ggplot(qc_highmaf,aes(x=missingness))+geom_histogram(position = "dodge",binwidth = 0.005,fill="#66C2A5") + labs(x="Proportion of missing data",y="Number of SNPs",title = "SNPs with MAF > 0.05")+
    geom_vline(xintercept=0.01, lty=2, col="red")+
    scale_color_brewer(palette = "Dark2")

  p_lowMAF <- ggplot(qc_lowmaf,aes(x=missingness))+geom_histogram(position = "dodge",binwidth=0.005,fill="#66C2A5") + labs(x="Proportion of missing data",y="Number of SNPs", title = "SNPs with MAF < 0.05") +
    geom_vline(xintercept=0.01, lty=2, col="red")+
    scale_color_brewer(palette = "Dark2")

  mafmiss <- cowplot::plot_grid(p_lowMAF, p_highMAF)


  #Minor allele frequency distribution
  maf <- ggplot(qc_report,aes(x=maf))+geom_histogram(binwidth= 0.01,fill="#66C2A5")+ labs(x="Minor allele frequency",y="Number of SNPs", title = "Minor allele frequency distribution")+ geom_vline(xintercept = 0.05, lty=2, col="red")

  #Hardy weinberg exact test p-val distribution
  qc_report$hwe_p_log <- -log10(qc_report$hwe_p)
  qc_lowhwe <- subset(qc_report,qc_report$hwe_p < 0.01)

  hwe_all <- ggplot(qc_report,aes(x=hwe_p_log))+geom_histogram(binwidth = 0.5,fill="#66C2A5")+ labs(x=expression("-log"[10]*" of HWE exact p-value"),y="Number of SNPs", title = "HWE exact")+ geom_vline(xintercept= 5, lty=2, col="red")
  hwe_low <- ggplot(qc_lowhwe,aes(x=hwe_p_log))+geom_histogram(binwidth = 0.5,fill="#66C2A5")+ labs(x=expression("-log"[10]* " of HWE exact p-value"),y="Number of SNPs", title = "HWE exact p-val <0.01")+ geom_vline(xintercept= 5, lty=2, col="red")

  hwes <- cowplot::plot_grid(hwe_all,hwe_low)

  plots_markerQC <- list(mafmiss, maf, hwes)
  subplotLabels <- LETTERS[1:length(plots_markerQC)]

  final_plot_all <- cowplot::plot_grid(plotlist=plots_markerQC,
                                   nrow=length(plots_markerQC),
                                   labels=subplotLabels,
                                   rel_heights=c(rep(1,
                                                     length(plots_markerQC)),
                                                 1.5))
  return(final_plot_all)
}


autoplot_l_qc_overview <- function(.x,...){


  qc_report <- .x

  qc_lowmaf <- subset(qc_report, qc_report$maf <= 0.05)
  qc_lowhwe <- subset(qc_report,qc_report$hwe_p < 0.01)


  maf_fail <- c(qc_lowmaf$snp_id)
  hwe_fail <- c(qc_lowhwe$snp_id)

  qc_highmissing <- subset(qc_report,qc_report$missingness>=0.05)
  missing_fail <- c(qc_highmissing$snp_id)

  list <- list(MAF = maf_fail, HWE = hwe_fail, Missing = missing_fail)

  final_plot_overview <- UpSetR::upset(UpSetR::fromList(list),order.by = "freq",main.bar.color="#66C2A5", matrix.color="#66C2A5",
                                       sets.bar.color="#FC8D62")

  return(final_plot_overview)
}
