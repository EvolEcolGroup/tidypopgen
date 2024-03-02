#' Filter, flip and swap SNPs in a genlight object
#'
#' This function takes a genlight object and edits its SNPs
#'
#' Note that this function does not check that to_keep and to_swap have the correct
#' length
#'
#' @param x a [adegenet::SNPbin] object.
#' @param new_order a numeric vector with the new order of the SNPs (NAs will
#' be excluded)
#' @param to_flip a boolean vector of SNPs to flip
#' @param to_swap a boolean vector of SNPs to swap
#' @param missing_allele a character vector of alternate alleles to fill in
#' missing alleles (NA for loci without missing alleles)
#' @returns an edited [adegenet::SNPbin] object.
#' @keywords internal
select_swap_snpbin_list <- function(x, new_order, to_flip, to_swap, missing_allele){
  # ordered ids without NAs
  new_order_na <- order(new_order,na.last=NA)
  # copy the loci attribute
  loci_info <- attr(x,"loci")
  # fix missing alleles
  loci_info$allele_alt[!is.na(missing_allele)]<-missing_allele[!is.na(missing_allele)]
  loci_info$allele_alt[to_flip] <- flip(loci_info$allele_alt[to_flip])
  loci_info$allele_ref[to_flip] <- flip(loci_info$allele_ref[to_flip])
  loci_info$allele_alt[to_swap] <- loci_info$allele_ref[to_swap]
  loci_info$allele_ref[to_swap] <- loci_info$allele_alt[to_swap]
  loci_info <- loci_info[new_order_na,]

  # now deal with the genotypes
  x <- lapply(x, select_swap_snpbin, sel_indices = new_order, to_swap = to_swap)
  attr(x,"loci") <- loci_info
  x
}
