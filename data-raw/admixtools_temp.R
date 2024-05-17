###########################################################################################
# everything works up to here
# the function below needs to be modified to use gt_extract_afs


gt_extract_f2_large = function(x, outdir,  blgsize = 0.05, cols_per_chunk = 10,
                               maxmiss = 0, minmaf = 0, maxmaf = 0.5, minac2 = FALSE, outpop = NULL, outpop_scale = TRUE,
                               transitions = TRUE, transversions = TRUE,
                               keepsnps = NULL, snpblocks = NULL, overwrite = FALSE, format = NULL,
                               adjust_pseudohaploid = TRUE, afprod = TRUE, fst = TRUE, poly_only = c('f2'),
                               apply_corr = TRUE, verbose = TRUE) {
  inds = NULL
  pops = NULL
  if(!quiet) message(paste0('Extracting allele frequencies...\n'))
  snpdat = extract_afs(pref, outdir, inds = inds, pops = pops, cols_per_chunk = cols_per_chunk, numparts = 100,
                       maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf, minac2 = minac2, outpop = outpop,
                       transitions = transitions, transversions = transversions,
                       keepsnps = keepsnps, format = format, poly_only = FALSE,
                       adjust_pseudohaploid = adjust_pseudohaploid, verbose = verbose)

  numchunks = length(list.files(outdir, 'afs.+rds'))

  if(!is.null(outpop) && outpop_scale) {
    p = snpdat$outpopaf
    snpwt = 1/(p*(1-p))
  } else snpwt = NULL
  if(isTRUE(poly_only)) poly_only = c('f2', 'ap', 'fst')

  if(verbose) alert_warning(paste0('Computing ', choose(numchunks+1, 2), ' chunk pairs. If this takes too long,
  consider running "extract_afs" and then paralellizing over "afs_to_f2".\n'))
  for(i in 1:numchunks) {
    for(j in i:numchunks) {
      if(verbose) alert_info(paste0('Writing pair ', i, ' - ', j, '...\r'))
      afs_to_f2(outdir, outdir, chunk1 = i, chunk2 = j, blgsize = blgsize, snpdat = snpdat,
                snpwt = snpwt, overwrite = overwrite, type = 'f2', poly_only = 'f2' %in% poly_only,
                apply_corr = apply_corr)
      if(afprod) afs_to_f2(outdir, outdir, chunk1 = i, chunk2 = j, blgsize = blgsize, snpdat = snpdat,
                           snpwt = snpwt, overwrite = overwrite, type = 'ap', poly_only = 'ap' %in% poly_only)
      if(fst) afs_to_f2(outdir, outdir, chunk1 = i, chunk2 = j, blgsize = blgsize, snpdat = snpdat,
                        snpwt = snpwt, overwrite = overwrite, type = 'fst', poly_only = 'fst' %in% poly_only)
    }
  }
  if(verbose) alert_info('\n')
  if(verbose) alert_info(paste0('Deleting allele frequency files...\n'))
  unlink(paste0(outdir, c('/afs*.rds', '/counts*.rds')))
}


## threre is something off with this function, ti does not seem to replicate extract_afs
gt_extract_afs <- function (.x, outdir = "./afs", cols_per_chunk = 10,
                            numparts = 100, maxmiss = 0, minmaf = 0, maxmaf = 0.5, minac2 = FALSE,
                            outpop = NULL, transitions = TRUE, transversions = TRUE,
                            format = NULL, poly_only = FALSE, adjust_pseudohaploid = TRUE,
                            quiet = FALSE){
  # variables that don't make sense with gen_tibble
  inds = NULL
  pops = NULL
  auto_only = FALSE
  keepsnps = NULL

  warning ("this function is untested")
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
  afdat %<>% discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf,
                                  maxmaf = maxmaf, minac2 = minac2, outpop = outpop, transitions = transitions,
                                  transversions = transversions, keepsnps = keepsnps, auto_only = auto_only,
                                  poly_only = poly_only)
  afdat$snpfile <- afdat$snpfile %>%  mutate(poly = as.logical(cpp_is_polymorphic(afdat$afs)))
  if (!quiet)
    message(paste0(nrow(afdat$afs), " SNPs remain after filtering. ",
                   sum(afdat$snpfile$poly), " are polymorphic.\n"))
  admixtools::split_mat(afdat$afs, cols_per_chunk = cols_per_chunk, prefix = paste0(outdir,
                                                                                    "/afs"), verbose = verbose)
  admixtools::split_mat(afdat$counts, cols_per_chunk = cols_per_chunk,
                        prefix = paste0(outdir, "/counts"), verbose = verbose)
  block_lengths = admixtools::get_block_lengths(afdat$snpfile %>% filter(.data$poly),
                                                blgsize = blgsize)
  block_lengths_a = admixtools::get_block_lengths(afdat$snpfile, blgsize = blgsize)
  saveRDS(block_lengths, file = paste0(outdir, "/block_lengths.rds"))
  saveRDS(block_lengths_a, file = paste0(outdir, "/block_lengths_a.rds"))
  readr::write_tsv(afdat$snpfile, paste0(outdir, "/snpdat.tsv.gz"))
  invisible(afdat$snpfile)
}

#' Compute and store blocked f2 statistics for a `gen_tibble`
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
#' @keywords internal

## TO BE TESTED
gt_extract_f2_broken <- function(.x , outdir=NULL, cols_per_chunk = 10, blgsize = 0.05, quiet=FALSE){
  warning ("this function is untested")
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
  gt_extract_afs(.x, outdir = af_dir, cols_per_chunk = cols_per_chunk,
                 blgsize = blgsize, quiet= quiet)
  numchunks = length(list.files(af_dir, 'afs.+rds'))
  for(i in 1:numchunks) {
    for(j in i:numchunks) {
      admixtools::afs_to_f2(af_dir, outdir, chunk1 = i, chunk2 = j)
    }
  }
}

#gt_extract_f2(test_gt)
