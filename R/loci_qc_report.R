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
  qc_report <- .x %>% reframe(maf=loci_freq(.data$genotypes),
                              missingness = loci_missingness(.data$genotypes),
                    hwe_p = loci_hwe(.data$genotypes, ...))
 class(qc_report) <- c("loci_qc_report",class(qc_report))
 qc_report
}

