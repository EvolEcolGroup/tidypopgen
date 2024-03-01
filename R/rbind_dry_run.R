#' Generate a report of what would happen to each SNP in a merge
#'
#' This function provides an overview of the fate of each SNP in two
#' [`adegenet::genlight`] objects in the case of a merge. Only SNPs found in
#' both objects will be kept. One object is used as a `reference`, and SNPs
#' in the other dataset will be flipped and/or alleles swapped
#' as needed. SNPs that have different alleles in the two datasets will also be
#' dropped.
#'
#' @param ref either a [`gen_tibble`] object, or the path to the PLINK file;
#' the alleles in this objects will
#' be used as template to flip the ones in `target` and/or
#' swap their order as necessary.
#' @param target either a [`gen_tibble`] object, or the path to the PLINK file
#' (saved in raw format, see details in [read_plink_raw()].
#' @param flip_strand boolean on whether strand flipping should be checked to
#' match the two datasets. It defaults to FALSE
#' @param remove_ambiguous boolean whether ambiguous SNPs (i.e. a/t and c/g)
#' should be removed. It defaults to FALSE
#' @param quiet boolean whether to omit reporting to screen
#' @returns a list with two `data.frames`, named `target` and `ref`. Each
#' data.frame has `nrow()` equal to the number of loci in the respective dataset,
#' a column `id` with the locus name, and boolean columns `to_keep` (the valid loci
#' that will be kept in the merge),
#' `alleles_mismatched` (loci found in both datasets but with mismatched alleles,
#' leading to those loci being dropped), `to_flip` (loci that need to be flipped
#' to align the two datasets, only found in `target` data.frame)
#' and `to_swap` (loci for which the order of alleles needs to be swapped
#' to align the two datasets, `target` data.frame)
#' @export

rbind_dry_run <- function(ref, target, flip_strand = FALSE,
                          remove_ambiguous = FALSE, quiet = FALSE){
  # create a data.frame with loci names, numeric_id, and alleles
  target_df <- target %>% show_loci()
  target_df <- target_df %>% mutate(id = seq_len(nrow(target_df)))
  ref_df <- ref %>% show_loci()
  ref_df <- ref_df %>% mutate(id = seq_len(nrow(ref_df)))
  # replace NA with "m" for missing allele to avoid subsetting headaches
  ref_df$allele_alt[is.na(ref_df$allele_alt)]<-"m"
  target_df$allele_alt[is.na(target_df$allele_alt)]<-"m"
  ref_df <- ref_df %>% rename(allele_1 = allele_alt, allele_2 = allele_ref)
  target_df <- target_df %>% rename(allele_1 = allele_alt, allele_2 = allele_ref)
  #browser()
  rbind_dry_run_df(ref_df = ref_df,
                   target_df = target_df,
                   flip_strand = flip_strand,
                     remove_ambiguous = remove_ambiguous,
                     quiet = quiet)
  }


##############################################################################
# Methods can simply call the function below after having correctly formatted
# the dataframes
##############################################################################

rbind_dry_run_df <- function(target_df, ref_df, flip_strand, remove_ambiguous, quiet){
  #browser()
  # now filter for alleles in common
  target_sub <- target_df[target_df$name %in% ref_df$name,]
  ref_sub <- ref_df[ref_df$name %in% target_df$name,]
  # reorder target_sub to match ref_sub
  target_sub <- target_sub[match(target_sub$name, ref_sub$name),]
  # we now have two data.frames with the same loci and in the same order

  # fix any missing alleles
  target_sub$missing_allele <- resolve_missing_alleles(missing_table = target_sub, other_table = ref_sub)
  target_missing_to_fix <- target_sub$allele_1=="m" & !is.na(target_sub$missing_allele)
  target_sub$allele_1[target_missing_to_fix] <- target_sub$missing_allele[target_missing_to_fix]
  ref_sub$missing_allele <- resolve_missing_alleles(missing_table = ref_sub, other_table = target_sub)
  ref_missing_to_fix <- ref_sub$allele_1=="m" & !is.na(ref_sub$missing_allele)
  ref_sub$allele_1[ref_missing_to_fix] <- ref_sub$missing_allele[ref_missing_to_fix]

  # preliminary list of alleles to keep based on whether alleles are correct
  # (irrespective of order, we will deal with that later)
  to_keep_orig <- (((target_sub$allele_1 == ref_sub$allele_1) &
                      (target_sub$allele_2 == ref_sub$allele_2)) |
                     (target_sub$allele_1 == ref_sub$allele_2) &
                     (target_sub$allele_2 == ref_sub$allele_1))

  if (flip_strand){
    # flip the ones that are not correct (in case it solves the problem)
    target_sub$allele_1[!to_keep_orig] <- flip(target_sub$allele_1[!to_keep_orig])
    target_sub$allele_2[!to_keep_orig] <- flip(target_sub$allele_2[!to_keep_orig])
  }


  # redefine alleles to keep (to include the ones we fixed with the flip)
  to_keep_flip <- (((target_sub$allele_1 == ref_sub$allele_1) &
                      (target_sub$allele_2 == ref_sub$allele_2)) |
                     (target_sub$allele_1 == ref_sub$allele_2) &
                     (target_sub$allele_2 == ref_sub$allele_1))

  to_flip <- (to_keep_flip & !to_keep_orig)
  # and now check which have to be swapped
  to_swap <- (target_sub$allele_1 == ref_sub$allele_2) &
    (target_sub$allele_2 == ref_sub$allele_1)
  if(remove_ambiguous){
    # remove ambiguous snps from the boolean vectors
    ambiguous_sub <- ambiguous(target_sub)
    to_keep_flip <- to_keep_flip & !ambiguous_sub
    to_flip <- to_flip & !ambiguous_sub
    to_swap <- to_swap & !ambiguous_sub
  }
  # now we create the two reporting data.frames
  # note that they include all loci (including the ones we dropped because they
  # did not exist in one of the datasets)
  # first create a report for the ref
  ref_report <- data.frame(id = ref_df$id,
                           new_id = match(ref_df$name, ref_sub$name[to_keep_flip]),
                           name = ref_df$name,
                           missing_allele = NA,
                           ambiguous = ambiguous (ref_df))
  ref_report$missing_allele[match(ref_sub$name,ref_df$name)] <- ref_sub$missing_allele
  # and now one for target
  target_report <- data.frame(id = target_df$id,
                              new_id = match(target_df$name, ref_sub$name[to_keep_flip]),
                              name = target_df$name,
                              to_flip = FALSE,
                              to_swap = FALSE,
                              missing_allele = NA,
                              ambiguous = ambiguous(target_df))

  # update the to_keep list
  target_report$to_flip[target_report$name %in% target_sub$name[to_flip]]<-TRUE
  target_report$to_swap[target_report$name %in% target_sub$name[to_swap]]<-TRUE
  target_report$missing_allele[match(target_sub$name,target_df$name)] <- target_sub$missing_allele

  report <- list(target = target_report, ref = ref_report)
  class(report) <- c("rbind_report",class(report))
  if (!quiet){
    summary(report)
  }
  return(report)
}



resolve_missing_alleles <- function(missing_table, other_table){
  #browser()
  missing_replacements <- rep(NA, nrow(missing_table))
  for (i_row in which(missing_table$allele_1=="m")){
    # get non-missing allele
    known_allele <- missing_table[i_row,"allele_2"]
    other_alleles <- unlist(other_table[i_row,c("allele_1","allele_2")])
    # as long as we don't have a missing allele in the other table as well
    if (other_alleles[1]!="m"){
      if (known_allele %in% other_alleles){
        missing_replacements[i_row] <- other_alleles[!other_alleles %in% known_allele]
      } else if (known_allele %in% flip(other_alleles)){
        missing_replacements[i_row] <- flip(other_alleles)[!flip(other_alleles) %in% known_allele]
      }
    }
    # if both datasets have missing alleles
    # we don't need to do anything
    # if (known_allele %in% other_allele | known_allele %in% flip(other_allele))
    # but what do we do if they don't match? Arguably just ditch the locus
  }
  return(missing_replacements)
}


ambiguous <- function(alleles_df){
  (((alleles_df$allele_1=="a") & (alleles_df$allele_2=="t")) |
     ((alleles_df$allele_1=="t") & (alleles_df$allele_2=="a")) |
     ((alleles_df$allele_1=="c") & (alleles_df$allele_2=="g")) |
     ((alleles_df$allele_1=="g") & (alleles_df$allele_2=="c")))
}

