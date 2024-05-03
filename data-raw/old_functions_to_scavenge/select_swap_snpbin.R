#' Select and swap SNPs in a SNP bin
#'
#' This function first swaps the major and minor alleles for loci set at TRUE
#' in to_swap, and then keeps
#'
#' Note that this function does not check that `to_swap` has the correct
#' length, nor that `sel_indices` contains valid indices (or it is a boolean
#' of the appropriate length). This is to optmise speed.
#'
#' @param x a [adegenet::SNPbin] object.
#' @param sel_indices the indices and order in which loci should be rearranged,
#' or a boolean of the appropriate length
#' @param to_swap indices or a boolean vector of SNPs to swap
#' @returns an edited [adegenet::SNPbin] object.
#' @keywords internal

select_swap_snpbin <- function(x, sel_indices, to_swap){
  # unpack the snp (can this be faster???)
  # TODO check the other SNPbin functions
  this_ploidy <- x@ploidy
  this_label <- x@label
  x <- as.integer(x)
  #.SNPbin2int(x)

  x[to_swap] <- dplyr::case_match(x[to_swap], 0 ~ 2, 2 ~ 0,
                                  .default = x[to_swap])
  methods::new("SNPbin", x[sel_indices], label=this_label, ploidy = this_ploidy)
}
