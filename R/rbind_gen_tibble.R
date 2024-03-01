#' @export
rbind.gen_tbl <- function(ref, target, flip_strand = FALSE,
              remove_ambiguous = FALSE, quiet = FALSE){
  report <- rbind_dry_run(target, ref, flip_strand=flip_strand,
                          remove_ambiguous = remove_ambiguous, quiet = quiet)
  # now edit the gen_tibble objects
  # for the ref object, we fix the missing alleles
  id_missing <- which(!is.na(report$ref$missing_allele))
  stop("this is not yet implemented!!!")
  if (length(id_missing)>0){
    adegenet::alleles(ref)[id_missing] <- paste0(report$ref$missing_allele[id_missing],
                                                 substr(adegenet::alleles(ref)[id_missing],2,3))
  }
  # and then use adegenet native cropping
  ref_to_keep <- which(!is.na(report$ref$new_id))
  ref_sub <- ref[,which(!is.na(report$ref$new_id))]
  # for the target, we use our own function as the transformation is more complicated
  target_sub <- filter_flip_genlight(target, new_order = report$target$new_id,
                                     to_flip = report$target$to_flip,
                                     to_swap = report$target$to_swap,
                                     missing_allele = report$target$missing_allele)
  other_table <- rbind(ref_sub@other$ind.metrics, target_sub@other$ind.metrics)
  ref_sub <- rbind(ref_sub,target_sub)
  adegenet::other(ref_sub)$ind.metrics <- other_table
  return(ref_sub)
}
