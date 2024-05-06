# admixtools equivalent functions

gt_to_aftable <- function(.x){
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  .group_levels = .x %>% group_keys() %>% pull(1)
  # summarise population frequencies
  pop_freqs_list <- group_map(.x, .f=~.gt_pop_freqs(.x))

  afs <- do.call("cbind",lapply(pop_freqs_list,function(x)x[["freq_alt"]]))
  dimnames(afs)<-list(loci_names(.x),.group_levels)
  counts <- do.call("cbind",lapply(pop_freqs_list,function(x)x[["n"]]))
  dimnames(counts)<-list(loci_names(.x),.group_levels)
  loci_new_names <- c(SNP = "name", CHR = "chromosome", POS="position",cm="genetic_dist",
                      A1 = "allele_ref",A2 = "allele_alt")
  snp <- show_loci(.x) %>% select(-all_of("big_index")) %>%
    rename(dplyr::all_of(loci_new_names)) %>% dplyr::relocate(all_of("cm"),.before="POS")
  c("SNP","CHR","cm","POS","A1","A2")

  aftable <- list(afs = afs,
                  counts = counts,
                  snpfile = snp)
  return(aftable)
}

gt_extract_afs <- function(.x, outdir = "./afs", cols_per_chunk = 10, blgsize = 0.05, quiet=FALSE){
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  if (!dir.exists(outdir)){
    dir.create(outdir)
  }
  verbose <- !quiet
  afdat = gt_to_aftable(.x)#, inds = inds, pops = pops,
                           #   format = format, adjust_pseudohaploid = adjust_pseudohaploid,
                           #   verbose = verbose)
#  afdat %<>% discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf,
#                                  maxmaf = maxmaf, minac2 = minac2, outpop = outpop, transitions = transitions,
#                                  transversions = transversions, keepsnps = keepsnps, auto_only = TRUE,
#                                  poly_only = poly_only)
  afdat$snpfile <- afdat$snpfile %>%  mutate(poly = as.logical(admixtools:::cpp_is_polymorphic(afdat$afs)))
#  if (verbose)
#    alert_warning(paste0(nrow(afdat$afs), " SNPs remain after filtering. ",
#                         sum(afdat$snpfile$poly), " are polymorphic.\n"))
  admixtools::split_mat(afdat$afs, cols_per_chunk = cols_per_chunk, prefix = paste0(outdir,
                                                                        "/afs"), verbose = verbose)
  admixtools::split_mat(afdat$counts, cols_per_chunk = cols_per_chunk,
            prefix = paste0(outdir, "/counts"), verbose = verbose)
  block_lengths = admixtools::get_block_lengths(afdat$snpfile %>% filter(poly),
                                    blgsize = blgsize)
  block_lengths_a = admixtools::get_block_lengths(afdat$snpfile, blgsize = blgsize)
  saveRDS(block_lengths, file = paste0(outdir, "/block_lengths.rds"))
  saveRDS(block_lengths_a, file = paste0(outdir, "/block_lengths_a.rds"))
  readr::write_tsv(afdat$snpfile, paste0(outdir, "/snpdat.tsv.gz"))
  invisible(afdat$snpfile)
}

#' Compute adn store blocked f2 statistics for a `gen_tibble`
#'
#' This function prepares data for various `ADMIXTOOLS 2` function, and it is
#' equivalent to `admixtools::extract_f2`. An important difference is that
#' the filtering for snps (e.g. by maf) or populations is not performed by this
#' function; it is expected that the `gen_tibble` has been filtered appropriately.
#'
#' @param .x the `gen_tibble`, appropriately filtered for individuals, populations and snps
#' @param outdir the directory where the f2 stats will be stored. If left NULL, a directory
#' names 'f2' will be created in the same path as the RDS fo the `gen_tibble`
#' @param blgsize SNP block size in Morgan. Default is 0.05 (5 cM). If `blgsize`
#' is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
#' @param cols_per_chunk Number of allele frequency chunks to store on disk.
#' Setting this to a positive integer makes the function slower,
#' but requires less memory. The default value for cols_per_chunk
#' in extract_afs is 10. Lower numbers will lower the memory requirement
#' but increase the time it takes.
#' @param quiet boolean on whether the progess should be silenced
#' @returns SNP metadata (invisibly)
#' @export

gt_extract_f2 <- function(.x , outdir=NULL, cols_per_chunk = 10, blgsize = 0.05, quiet=FALSE){
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  # if no outdir is given, create a subdirectory f2 in the path of the gen_tibble rds
  if (is.null(outdir)){
    outdir <- file.path(dirname(.gt_get_bigsnp(.x)$genotypes$rds),"f2")
  }
  # we will place the af tables into a subdirectory
  af_dir <- file.path(outdir,"af_tbls")
  if (!dir.exists(af_dir)){
    dir.create(af_dir,recursive = TRUE)
  }
  gt_extract_afs(test_gt, outdir = af_dir, cols_per_chunk = cols_per_chunk,
                 blgsize = blgsize, quiet= quiet)
  numchunks = length(list.files(af_dir, 'afs.+rds'))
  for(i in 1:numchunks) {
    for(j in i:numchunks) {
      admixtools::afs_to_f2(af_dir, outdir, chunk1 = i, chunk2 = j)
    }
  }
}

#gt_extract_f2(test_gt)
