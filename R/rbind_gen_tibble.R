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
#' @param backingfile the path and prefix of the files used to store the
#' merged data (it will be a .RDS to store the `bigSNP` object and a .bk file
#' as its backing file for the FBM)
#' @returns a [`gen_tibble`] with the merged data.
#' @export
rbind.gen_tbl <- function(..., as_is = FALSE, flip_strand = FALSE,
              remove_ambiguous = FALSE, quiet = FALSE, backingfile=NULL){
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
  # sort out paths
  if (is.null(backingfile)){
    save_path <- dirname(attr(ref$genotypes,"bigsnp")$genotypes$backingfile)
    backingfile <- tempfile("gt_merged_",tmpdir = save_path, fileext = "")
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
  # and in  the bigSNP object
  attr(ref$genotypes,"bigsnp")$map$allele1[attr(ref$genotypes, "loci")$big_index[id_missing]] <- report$ref$missing_allele[id_missing]
  # now create a new loci table (we'll use it later)
  new_ref_loci_tbl <- show_loci(ref)[order(report$ref$new_id,na.last=NA),]
  # now we subset the SNP object
  ## in the snp object
  ref_snp <- subset_bigSNP(attr(ref$genotypes,"bigsnp"),
                           loci_indices = new_ref_loci_tbl$big_index,
                           indiv_indices = vctrs::vec_data(ref$genotypes))
  ###########
  # now we move to the target object
  # we fix the missing alleles
  ## in the gt_table
  id_missing <- which(!is.na(report$target$missing_allele))
  attr(target$genotypes, "loci")$allele_alt[id_missing] <- report$target$missing_allele[id_missing]
  # and in  the bigSNP object
  attr(target$genotypes,"bigsnp")$map$allele1[attr(target$genotypes, "loci")$big_index[id_missing]] <- report$target$missing_allele[id_missing]
  # now flip the strands
  ## in the gt_table

  to_flip <- report$target$to_flip
  attr(target$genotypes, "loci")$allele_alt[to_flip] <- flip(attr(target$genotypes, "loci")$allele_alt[to_flip])
  attr(target$genotypes, "loci")$allele_ref[to_flip] <- flip(attr(target$genotypes, "loci")$allele_ref[to_flip])
  to_flip_big_index <- show_loci(target)$big_index[to_flip]
  attr(target$genotypes,"bigsnp")$map$allele1[to_flip_big_index] <- flip(attr(target$genotypes,"bigsnp")$map$allele1[to_flip_big_index])
  attr(target$genotypes,"bigsnp")$map$allele2[to_flip_big_index] <- flip(attr(target$genotypes,"bigsnp")$map$allele2[to_flip_big_index])

  # now create a new loci table (we'll use it later)
  new_target_loci_tbl <- show_loci(target)[order(report$target$new_id,na.last=NA),]
  # now we subset the SNP object
  ## in the snp object
  target_snp <- subset_bigSNP(attr(target$genotypes,"bigsnp"),
                           loci_indices = new_target_loci_tbl$big_index,
                           indiv_indices = vctrs::vec_data(target$genotypes),
                           swap_indices = show_loci(target)$big_index[report$target$to_swap])

  if (!identical(target_snp$map$marker.ID,ref_snp$map$marker.ID) |
      !identical(target_snp$map$allele2,ref_snp$map$allele2)) {
    stop("something went wrong when subsetting and reordering the two bigSNP objects")
  }

  # now we need to merge the two FBMs
  # we start by transposing them, so that they just need to be appended to each other
  t_ref_fbm <- bigstatsr::big_transpose(ref_snp$genotypes)
  t_target_fbm <- bigstatsr::big_transpose(target_snp$genotypes)
  # append the two files
  append_success <- file.append(t_ref_fbm$backingfile, t_target_fbm$backingfile)
  # and amend the new number of columns
  t_ref_fbm$ncol<-t_ref_fbm$ncol+t_target_fbm$ncol
  # now flip the file around
  merged_fbm <- bigstatsr::big_transpose(t_ref_fbm, backingfile = backingfile) # TODO this should be written in the director of interest
  merged_snp <- structure(list(genotypes = merged_fbm,
                             fam = rbind(ref_snp$fam, target_snp$fam),
                             map = ref_snp$map),
                        class = "bigSNP")
  merged_rds <- paste0(backingfile,".rds")
  saveRDS(merged_snp, merged_rds)
  # Now we need to create the gen_tibble
  # TODO we need to turn genotype into a string!!!
  merged_tbl <- rbind(ref %>% select(-dplyr::any_of("genotypes")), target %>% select(-dplyr::any_of("genotypes")))
  # make sure that the genotypes vector points to the correct rows
  vctrs::vec_data(ref$genotypes)
  #and finally append the loci table
  indivs_with_big_names <- c(names(ref$genotypes),names(target$genotypes))
  #browser()
  new_ref_loci_tbl$big_index<-match(new_ref_loci_tbl$name,merged_snp$map$marker.ID) # TODO check that this is the correct order!!!!
  # TODO check that all individuals in tibble and bigsnp object are the same
  merged_tbl$genotypes <- vctrs::new_vctr(match(indivs_with_big_names,merged_snp$fam$sample.ID), # TODO check that this is the correct order!!!!
                  bigsnp = merged_snp,
                  loci=new_ref_loci_tbl,
                  names=indivs_with_big_names,
                  class = "vctrs_bigSNP")

  # TODO check that the snp is saved (it should be), and let the user know what the new
  # name of the files are!!!
  if (!quiet){
    message("\n\nthe new bigSNP file is: ", merged_rds)
    message("the new backing file is: ", attr(merged_tbl$genotypes,"bigsnp")$genotypes$backingfile)
  }
  tibble::new_tibble(
    merged_tbl,
    class = "gen_tbl"
  )
}

