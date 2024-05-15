#' Compute and store blocked f2 statistics
#'
#' This function prepares data for various other *ADMIXTOOLS 2* functions. It
#' takes a [`gen_tibble`],
#' computes allele frequencies and blocked f2-statistics for selected populations,
#' and writes the results to `outdir`.
#' @export
#' @param .x a [`gen_tibble`]
#' @param outdir Directory where data will be stored.
#' @param blgsize SNP block size in Morgan. Default is 0.05 (5 cM). If `blgsize` is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
#' @param maxmem Maximum amount of memory to be used. If the required amount of memory exceeds `maxmem`, allele frequency data will be split into blocks, and the computation will be performed separately on each block pair. This doesn't put a precise cap on the amount of memory used (it used to at some point). Set this parameter to lower values if you run out of memory while running this function. Set it to higher values if this function is too slow and you have lots of memory.
#' @param maxmiss Discard SNPs which are missing in a fraction of populations higher than `maxmiss`
#' @param minmaf Discard SNPs with minor allele frequency less than `minmaf`
#' @param maxmaf Discard SNPs with minor allele frequency greater than than `maxmaf`
#' @param minac2 Discard SNPs with allele count lower than 2 in any population (default `FALSE`). This option should be set to `TRUE` when computing f3-statistics where one population consists mostly of pseudohaploid samples. Otherwise heterozygosity estimates and thus f3-estimates can be biased. `minac2 == 2` will discard SNPs with allele count lower than 2 in any non-singleton population (this option is experimental and is based on the hypothesis that using SNPs with allele count lower than 2 only leads to biases in non-singleton populations). While the `minac2` option discards SNPs with allele count lower than 2 in any population, the \code{\link{qp3pop}} function will only discard SNPs with allele count lower than 2 in the first (target) population (when the first argument is the prefix of a genotype file).
#' @param outpop Keep only SNPs which are heterozygous in this population
#' @param outpop_scale Scale f2-statistics by the inverse `outpop` heteroygosity (`1/(p*(1-p))`). Providing `outpop` and setting `outpop_scale` to `TRUE` will give the same results as the original *qpGraph* when the `outpop` parameter has been set, but it has the disadvantage of treating one population different from the others. This may limit the use of these f2-statistics for other models.
#' @param transitions Set this to `FALSE` to exclude transition SNPs
#' @param transversions Set this to `FALSE` to exclude transversion SNPs
#' @param auto_only Keep only SNPs on chromosomes 1 to 22
#' @param keepsnps SNP IDs of SNPs to keep. Overrides other SNP filtering options
#' @param overwrite Overwrite existing files in `outdir`
#' @param adjust_pseudohaploid Genotypes of pseudohaploid samples are usually coded as `0` or `2`, even though only one allele is observed. `adjust_pseudohaploid` ensures that the observed allele count increases only by `1` for each pseudohaploid sample. If `TRUE` (default), samples that don't have any genotypes coded as `1` among the first 1000 SNPs are automatically identified as pseudohaploid. This leads to slightly more accurate estimates of f-statistics. Setting this parameter to `FALSE` treats all samples as diploid and is equivalent to the *ADMIXTOOLS* `inbreed: NO` option. Setting `adjust_pseudohaploid` to an integer `n` will check the first `n` SNPs instead of the first 1000 SNPs.
#' @param cols_per_chunk Number of allele frequency chunks to store on disk. Setting this to a positive integer makes the function slower, but requires less memory. The default value for `cols_per_chunk` in \code{\link{extract_afs}} is 10. Lower numbers will lower the memory requirement but increase the time it takes.
#' @param fst Write files with pairwise FST for every population pair. Setting this to FALSE can make `extract_f2` faster and will require less memory.
#' @param afprod  Write files with allele frequency products for every population pair. Setting this to FALSE can make `extract_f2` faster and will require less memory.
#' @param poly_only Specify whether SNPs with identical allele frequencies in every population should be discarded (`poly_only = TRUE`), or whether they should be used (`poly_only = FALSE`). By default (`poly_only = c("f2")`), these SNPs will be used to compute FST and allele frequency products, but not to compute f2 (this is the default option in the original ADMIXTOOLS).
#' @param apply_corr Apply small-sample-size correction when computing f2-statistics (default `TRUE`)
#' @param n_cores Parallelize computation across `n_cores` cores via the `doParallel` package.
#' @param quiet Suppress printin of progress updates
#' @return SNP metadata (invisibly)
#' @export

gt_extract_f2<- function(.x , outdir=NULL, blgsize = 0.05, maxmem = 8000,
                         maxmiss = 0, minmaf = 0, maxmaf = 0.5, minac2 = FALSE, outpop = NULL, outpop_scale = TRUE,
                         transitions = TRUE, transversions = TRUE,
                         auto_only = TRUE, keepsnps = NULL, overwrite = FALSE, format = NULL,
                         adjust_pseudohaploid = TRUE, cols_per_chunk = NULL, fst = TRUE, afprod = TRUE,
                         poly_only = c('f2'), apply_corr = TRUE, qpfstats = FALSE, n_cores = 1, quiet = FALSE) {

  inds = NULL
  pops = NULL
  pops2 = NULL
  if (!is.null(cols_per_chunk)){
    stop("there is no option at the moment to split the processing in chunks")
  }
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  # if no outdir is given, create a subdirectory f2 in the path of the gen_tibble rds
  if (is.null(outdir)){
    outdir <- file.path(dirname(.gt_get_bigsnp(.x)$genotypes$rds),"f2")
  }

  verbose <- !quiet
  afdat = gt_to_aftable(.x)


  if(is.null(inds)) pops = union(pops, pops2)
  afdat %<>% admixtools:::discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf, minac2 = minac2, outpop = outpop,
                                               transitions = transitions, transversions = transversions,
                                               keepsnps = keepsnps, auto_only = auto_only, poly_only = FALSE)
  afdat$snpfile %<>% mutate(poly = as.logical(admixtools:::cpp_is_polymorphic(afdat$afs)))
  if(sum(afdat$snpfile$poly) == 0) stop('There are no informative SNPs!')

  if(verbose) message(paste0(nrow(afdat$afs), ' SNPs remain after filtering. ',
                             sum(afdat$snpfile$poly),' are polymorphic.\n'))

  if(isTRUE(poly_only)) poly_only = c('f2', 'ap', 'fst')
  arrs = admixtools:::afs_to_f2_blocks(afdat, outdir = outdir, overwrite = overwrite, maxmem = maxmem, poly_only = poly_only,
                                       pops1 = pops, pops2 = pops2, outpop = if(outpop_scale) outpop else NULL,
                                       blgsize = blgsize, afprod = afprod, fst = fst, apply_corr = apply_corr,
                                       n_cores = n_cores, verbose = verbose)

  if(is.null(outdir)) return(arrs)

  if(verbose) message(paste0('Data written to ', outdir, '/\n'))
  invisible(afdat$snpfile)
}



# admixtools equivalent functions
gt_to_aftable <- function(.x, n_cores = bigstatsr::nb_cores()){
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  geno_fbm <- .gt_get_bigsnp(.x)$genotypes

  aftable <- gt_grouped_alt_freq_diploid(BM = geno_fbm,rowInd = .gt_bigsnp_rows(.x),
                                          colInd = .gt_bigsnp_cols(.x),
                                          groupIds = dplyr::group_indices(.x)-1,
                                          ngroups = max(dplyr::group_indices(.x)),
                                          ncores = n_cores)
  names(aftable)<-c("afs","counts")

  .group_levels = .x %>% group_keys() %>% pull(1)
  dimnames(aftable$afs)<-list(loci_names(.x),.group_levels)
  dimnames(aftable$counts)<-list(loci_names(.x),.group_levels)

  loci_new_names <- c(SNP = "name", CHR = "chromosome", POS="position",cm="genetic_dist",
                      A1 = "allele_ref",A2 = "allele_alt")
  snp <- show_loci(.x) %>% select(-all_of("big_index")) %>%
    rename(dplyr::all_of(loci_new_names)) %>% dplyr::relocate(all_of("cm"),.before="POS")
  c("SNP","CHR","cm","POS","A1","A2")

  aftable$snpfile <- snp
  return(aftable)
}

###########################################################################################
# everything works up to here
# the function below needs to be modified to use gt_extract_afs


gt_extract_f2_large = function(pref, outdir, inds = NULL, pops = NULL, blgsize = 0.05, cols_per_chunk = 10,
                            maxmiss = 0, minmaf = 0, maxmaf = 0.5, minac2 = FALSE, outpop = NULL, outpop_scale = TRUE,
                            transitions = TRUE, transversions = TRUE,
                            keepsnps = NULL, snpblocks = NULL, overwrite = FALSE, format = NULL,
                            adjust_pseudohaploid = TRUE, afprod = TRUE, fst = TRUE, poly_only = c('f2'),
                            apply_corr = TRUE, verbose = TRUE) {

  if(verbose) alert_info(paste0('Extracting allele frequencies...\n'))
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
















gt_extract_afs <- function(.x, outdir = "./afs", cols_per_chunk = 10, blgsize = 0.05, quiet=FALSE){
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



