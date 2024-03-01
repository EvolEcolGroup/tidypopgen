#' @export
rbind_gen_tbl <- function(ref, target, flip_strand = FALSE,
              remove_ambiguous = FALSE, quiet = FALSE){
  report <- rbind_dry_run(target, ref, flip_strand=flip_strand,
                          remove_ambiguous = remove_ambiguous, quiet = quiet)
  # now edit the gen_tibble objects
  # for the ref object, we fix the missing alleles
  id_missing <- which(!is.na(report$ref$missing_allele))
  if (length(id_missing)>0){
    attr(ref, "loci")$allele_alt <- report$ref$missing_allele[id_missing]
#    adegenet::alleles(ref)[id_missing] <- paste0(report$ref$missing_allele[id_missing],
#                                                 substr(adegenet::alleles(ref)[id_missing],2,3))
  }
 # stop("this is not yet implemented!!!")

  # and then use adegenet native cropping
  # @TODO get the native cropping for the list of snpbin
#  ref_to_keep <- which(!is.na(report$ref$new_id))
#  ref_sub <- ref[,which(!is.na(report$ref$new_id))]
  ref$genotypes <- select_swap_snpbin_list(ref$genotypes, new_order = report$ref$new_id,
                                          to_flip = rep(FALSE, nrow(report$ref)),
                                          to_swap = rep(FALSE, nrow(report$ref)),
                                          missing_allele = rep(NA, nrow(report$ref)))
  # for the target, we use our own function as the transformation is more complicated
  target$genotypes <- select_swap_snpbin_list(target$genotypes, new_order = report$target$new_id,
                                     to_flip = report$target$to_flip,
                                     to_swap = report$target$to_swap,
                                     missing_allele = report$target$missing_allele)
  ref_sub <- rbind(ref,target)
  return(ref_sub)
}
