#' Clump loci based on a Linkage Disequilibrium threshold
#'
#' This function uses clumping to remove SNPs at high LD. When used with its
#' default options, clumping based on MAF is similar to standard pruning (as
#' done by PLINK with "--indep-pairwise (size+1) 1 thr.r2", but it results in a
#' better spread of SNPs over the chromosome.
#'
#' Any missing values in the genotypes of a `gen_tibble` passed to
#' `loci_ld_clump` will cause an error. To deal with missingness, see
#' [gt_impute_simple()].
#'
#' @param .x a [`gen_tibble`] object
#' @param thr_r2 Threshold over the squared correlation between two SNPs.
#'   Default is `0.2`.
#' @param S A vector of loci statistics which express the importance of each SNP
#'   (the more important is the SNP, the greater should be the corresponding
#'   statistic).\cr For example, if `S` follows the standard normal
#'   distribution, and "important" means significantly different from 0, you
#'   must use `abs(S)` instead.\cr
#' **If not specified, MAFs are computed and used.**
#' @param size For one SNP, window size around this SNP to compute correlations.
#'   Default is `100 / thr.r2` for clumping (0.2 -> 500; 0.1 -> 1000; 0.5 ->
#'   200). If `use_positions = FALSE`, this is a window in number of SNPs,
#'   otherwise it is a window in kb (genetic distance). Ideally, use positions,
#'   as they provide a more sensible approach.
#' @param use_positions boolean, if TRUE (the default), `size` is in kb, if
#'   FALSE size is the number of SNPs.
#' @param exclude Vector of SNP indices to exclude anyway. For example, can be
#'   used to exclude long-range LD regions (see Price2008). Another use can be
#'   for thresholding with respect to p-values associated with `S`.
#' @param n_cores number of cores to be used
#' @param return_id boolean on whether the id of SNPs to keep should be
#'   returned. It defaults to FALSE, which returns a vector of booleans (TRUE or
#'   FALSE)
#' @param ... currently not used.
#' @return a boolean vector indicating whether the SNP should be kept (if
#'   'return_id = FALSE', the default), else a vector of SNP indices to be kept
#'   (if 'return_id = TRUE')
#' @rdname loci_ld_clump
#' @export
loci_ld_clump <- function(.x, ...) {
  UseMethod("loci_ld_clump", .x)
}

#' @export
#' @rdname loci_ld_clump
loci_ld_clump.tbl_df <- function(.x, ...) {
  # TODO this is a hack to deal with the class being dropped when going through
  # group_map
  stopifnot_gen_tibble(.x)
  # check n_cores are available

  loci_ld_clump(.x$genotypes, ...)
}


#' @export
#' @rdname loci_ld_clump
loci_ld_clump.vctrs_bigSNP <- function(
    .x,
    S = NULL, # nolint
    thr_r2 = 0.2,
    size = 100 / thr_r2,
    exclude = NULL,
    use_positions = TRUE,
    n_cores = 1,
    return_id = FALSE,
    ...) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)

  if (n_cores > 1) {
    # Remove checking for two levels of parallelism
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = TRUE), add = TRUE)
  }

  # check that the loci have not been resorted
  # check that big_index in the loci table is an increasing sequence of indeces
  if (is.unsorted(show_loci(.x)$big_index, strictly = TRUE)) {
    stop(paste(
      "Your loci have been resorted; first save the new file backed",
      "matrix with `gt_update_baking_file()"
    ))
  }

  if (gt_has_imputed(.x) && gt_uses_imputed(.x) == FALSE) {
    # not uses_imputed
    gt_set_imputed(.x, set = TRUE)
    on.exit(gt_set_imputed(.x, set = FALSE), add = TRUE)
  }

  is_loci_table_ordered(.x, error_on_false = TRUE)

  # get the FBM
  geno_fbm <- attr(.x, "bigsnp")$genotypes # nolint
  # rows (individuals) that we want to use
  if (use_positions) {
    .positions <- rep(NA, nrow(attr(.x, "bigsnp")$map))
    .positions[show_loci(.x)$big_index] <- show_loci(.x)$position
  } else {
    .positions <- NULL
  }
  # create a chromosome vector (fill gaps between bigsnpr and show_loci)
  .chromosome <- rep(2147483647L, nrow(attr(.x, "bigsnp")$map))
  .chromosome[show_loci(.x)$big_index] <- show_loci(.x)$chr_int
  # now figure out if we have any snp which have already been removed
  # those will go into `exclude`
  loci_not_in_tibble <-
    seq_len(nrow(attr(.x, "bigsnp")$map))[
      !seq_len(nrow(attr(.x, "bigsnp")$map)) %in% .gt_bigsnp_cols(.x)
    ] # nolint
  exclude <- c(loci_not_in_tibble, .gt_bigsnp_cols(.x)[exclude])
  if (length(exclude) == 0) {
    exclude <- NULL
  }

  # as long as we have more than one individual
  snp_clump_ids <- bigsnpr::snp_clumping(
    G = attr(.x, "bigsnp")$genotypes,
    # infos.chr = show_loci(.x)$chr_int, #nolint start
    # TEMP HACK using the info from the bigsnpr object
    # infos.chr = cast_chromosome_to_int(attr(.x,"bigsnp")$map$chromosome), #nolint end
    infos.chr = .chromosome,
    ind.row = vctrs::vec_data(.x),
    S = S, # nolint
    thr.r2 = thr_r2,
    infos.pos = .positions,
    size = size,
    exclude = exclude,
    ncores = n_cores
  )
  to_keep_id <- match(snp_clump_ids, show_loci(.x)$big_index)

  if (return_id) {
    to_keep_id
  } else {
    keep_bool <- rep(FALSE, count_loci(.x))
    keep_bool[to_keep_id] <- TRUE
    keep_bool
  }
}

#' @export
#' @rdname loci_ld_clump
loci_ld_clump.grouped_df <- function(.x, ...) {
  # TODO this is seriously inefficient, we need to cast it into a
  # big_apply problem
  # of maybe it isn't that bad...
  # group_map(.x, .f=~loci_ld_clump(.x, ...)) #nolint start

  # ungroup .x #nolint end
  .x <- ungroup(.x)

  # pass to loci_ld_clump.tbl_df
  loci_ld_clump(.x, ...)
}
