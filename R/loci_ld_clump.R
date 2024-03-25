#' Clump loci based on a Linkage Disequilibrium threshold
#'
#' This function uses clumping to remove SNPs at high LD. When used with its
#' default options, clumping based on MAF is similar to standard pruning (as done by
#' PLINK with "--indep-pairwise (size+1) 1 thr.r2", but it results in a better
#' spread of SNPs over the chromosome.
#'
#' TODO we should really return a boolean rather than indices, so that it can
#' be easily used with `select_loci_if`
#'
#' @param .x a [`gen_tibble`] object
#' @param thr_r2 Threshold over the squared correlation between two SNPs.
#'   Default is `0.2`.
#' @param S A vector of loci statistics which express the importance
#' of each SNP (the more important is the SNP, the greater should be
#' the corresponding statistic).\cr
#' For example, if `S` follows the standard normal distribution, and "important"
#' means significantly different from 0, you must use `abs(S)` instead.\cr
#' **If not specified, MAFs are computed and used.**
#' @param size For one SNP, window size around this SNP to compute correlations.
#' Default is `100 / thr.r2` for clumping (0.2 -> 500; 0.1 -> 1000; 0.5 -> 200).
#' If `use_positions = FALSE`, this is a window in
#' number of SNPs, otherwise it is a window in kb (genetic distance).
#' Ideally, use positions, as they provide a more sensible approach.
#' @param use_positions boolean, if TRUE (the default), `size` is in kb, if FALSE
#' size is the number of SNPs.
#' @param exclude Vector of SNP indices to exclude anyway. For example,
#' can be used to exclude long-range LD regions (see Price2008). Another use
#' can be for thresholding with respect to p-values associated with `S`.
#' @param n_cores number of cores to be used
#' @param ... currently not used.
#' @return a vector of snp indices to be kept
#' @rdname loci_ld_clump
#' @export
loci_ld_clump <- function(.x, ...) {
    UseMethod("loci_ld_clump", .x)
}

#' @export
#' @rdname loci_ld_clump
loci_ld_clump.tbl_df <- function(.x, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_ld_clump(.x$genotypes, ...)
}


#' @export
#' @rdname loci_ld_clump
loci_ld_clump.vctrs_bigSNP <- function(.x,
                                       S = NULL,
                                       thr_r2 = 0.2,
                                       size = 100/thr_r2,
                                       exclude = NULL,
                                       use_positions = TRUE,
                                       n_cores = 1,
                                       ...)
{
  rlang::check_dots_empty()
  # get the FBM
  geno_fbm <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  if (use_positions){
    .positions <- show_loci(.x)$position
  } else {
    .positions <- NULL
  }
  # as long as we have more than one individual
  snp_clump_ids <- bigsnpr::snp_clumping(G = .x,
                        infos.chr = show_loci(.x)$chromosome,
                        ind.row = .gt_bigsnp_rows(.x),
                        ind.col = .gt_bigsnp_cols(.x),
                        S = S,
                        thr.r2 = thr_r2,
                        infos.pos = .positions,
                        size = size,
                        exclude = .gt_bigsnp_cols(.x)[exclude],
                        ncores = n_cores)
  warning("this is yet to be tested!!!")
  match(snp_clump_ids, show_loci(.x)$bid_id)
  ## @TODO test that this works as expected
}

#' @export
#' @rdname loci_ld_clump
loci_ld_clump.grouped_df <- function(.x, ...) {
  # TODO this is seriously inefficient, we need to cast it into a big_apply problem
  # of maybe it isn't that bad...
  group_map(.x, .f=~loci_ld_clump(.x, ...))
}



