#' Combine two gen_tibbles
#'
#' This function combined two [gen_tibble]s. By defaults, it subsets the loci
#' and swaps ref and alt alleles to make the two datasets compatible (this
#' behaviour can be switched off with `as_is`).
#' The first object is used as a "reference" , and SNPs
#' in the other dataset will be flipped and/or alleles swapped
#' as needed. SNPs that have different alleles in the two datasets
#' (i.e. triallelic) will also be
#' dropped. There are also options (NOT default) to attempt strand flipping to
#' match alleles (often needed in human datasets from different SNP chips),
#' and remove ambiguous alleles (c/g and a/t) where the correct strand can not
#' be guessed.
#'
#' @param ... two [`gen_tibble`] objects. Note tha this function can not take
#' more objects, `rbind` has to be done sequentially for large sets of objects.
#' @param as_is boolean determining whether the loci should be left as they are
#' before merging. If FALSE (the defaults), `rbind` will attempt to subset and
#' swap alleles as needed.
#' @param flip_strand boolean on whether strand flipping should be checked to
#' match the two datasets. It defaults to FALSE
#' @param remove_ambiguous boolean whether ambiguous SNPs (i.e. A/T and C/G)
#' should be removed. It defaults to FALSE
#' @param quiet boolean whether to omit reporting to screen
#' @returns a [`gen_tibble`] with the merged data.
#' @export
rbind_gen_tbl <- function(..., as_is = FALSE, flip_strand = FALSE,
              remove_ambiguous = FALSE, quiet = FALSE){
  dots <- list(...)
  if (length(dots)!=2){
    stop("rbind for gen_tibble can only take two tibbles at a time")
  }
  ref <-dots[[1]]
  target <- dots[[2]]
  if (!quiet){
    if (as_is){
      if (any(flip_strand, remove_ambiguous)){
        stop("if 'as_is' is set to TRUE, 'flip_strand' and 'remove_ambiguous' have to be FALSE")
      }
      message("rbind will not attempt to harmonise the loci (e.g. flip strand, reorder or subset\n",
              "set 'as_is' to FALSE to subset to loci in common\n")
    }
  }
  report <- rbind_dry_run(ref = ref, target = target, flip_strand=flip_strand,
                          remove_ambiguous = remove_ambiguous, quiet = quiet)
  # now edit the gen_tibble objects
  ###########
  # we start with the ref object
  # we fix the missing alleles
  ## in the gt_table
  id_missing <- which(!is.na(report$ref$missing_allele))
  attr(ref$genotypes, "loci")$allele_alt[id_missing] <- report$ref$missing_allele[id_missing]
  ## and then in the bigSNP object
  attr(ref$genotypes,"bigsnp")$map$allele
  # now we subset the SNP object
  ## in the gt_table

  ref <- ref %>% select_loci(order(report$ref$new_id,na.last=NA))

  # Now subset the FBM with a deep copy with indices

  # update the other tables in the bigSNP object

  # update the loci table in the gen_tibble

  ###########
  # now we move to the target object



  #
  # fix missing alleles in target
  id_missing_target <- which(!is.na(report$target$missing_allele))
  attr(target$genotypes, "loci")$allele_alt[id_missing_target] <- report$target$missing_allele[id_missing_target]
  to_flip <- report$target$to_flip
  attr(target$genotypes, "loci")$allele_alt[to_flip] <- flip(attr(target$genotypes, "loci")$allele_alt[to_flip])
  attr(target$genotypes, "loci")$allele_ref[to_flip] <- flip(attr(target$genotypes, "loci")$allele_ref[to_flip])

  ## TODO we should also flip to later check integrity before the merge
  target <- target %>% select_loci (order(report$target$new_id,na.last=NA),
                                    .swap_if_arg = report$target$to_swap)

  #check that the loci tables are the same
  if(!identical(show_loci(ref)%>% dplyr::select("name", "allele_ref", "allele_alt"),
            show_loci(target) %>% dplyr::select("name", "allele_ref", "allele_alt"))){
    if (as_is){
      stop("the two datasets have different loci; set 'as_is'=FALSE to let 'rbind' attempt to harmonise them")
    } else {
      stop("'rbind' attempted to harmonise the loci but something went WRONG")
    }
  }

  # for the use of the data.frame method (we can't use NextMethod with rbind)
  return(base::rbind.data.frame(ref,target))
}


subset_gen_tbl <- function(x, indiv_indices=NULL, loci_indices=NULL, swap_indices=NULL,
         backingfile = NULL) {

}
