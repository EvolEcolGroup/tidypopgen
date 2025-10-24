#' Clump loci based on a Linkage Disequilibrium threshold
#'
#' This function uses clumping to remove SNPs at high LD. When used with its
#' default options, clumping based on MAF is similar to standard pruning (as
#' done by PLINK with "--indep-pairwise (size+1) 1 thr.r2", but it results in a
#' better spread of SNPs over the chromosome. This function is a wrapper around
#' [bigsnpr::snp_clumping()]. See
#' https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html for more
#' information on the differences between pruning and clumping.
#'
#' Any missing values in the genotypes of a `gen_tibble` passed to
#' `loci_ld_clump` will cause an error. To deal with missingness, see
#' [gt_impute_simple()].
#'
#' @param .x a [`gen_tibble`] object
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#' directly to the function). This defaults to "genotypes" and can only take
#' that value. There is no need for the user to set it, but it is included to
#' resolve certain tidyselect operations.
#' @param thr_r2 Threshold over the squared correlation between two SNPs.
#'   Default is `0.2`.
#' @param S A vector of loci statistics which express the importance of each SNP
#'   (the more important is the SNP, the greater should be the corresponding
#'   statistic).\cr For example, if `S` follows the standard normal
#'   distribution, and "important" means significantly different from 0, you
#'   must use `abs(S)` instead.\cr
#' **If not specified, MAFs are computed and used.**
#' @param size For one SNP, window size around this SNP to compute correlations.
#'   Default is `100 / thr_r2` for clumping (0.2 -> 500; 0.1 -> 1000; 0.5 ->
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
#' @seealso [bigsnpr::snp_clumping()] which this function wraps.
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("gen_tbl") %>% gt_impute_simple()
#'
#' # To return a boolean vector indicating whether the SNP should be kept
#' example_gt %>% loci_ld_clump()
#' # To return a vector of SNP indices to be kept
#' example_gt %>% loci_ld_clump(return_id = TRUE)
#'
loci_ld_clump <- function(.x, .col = "genotypes", ...) {
  UseMethod("loci_ld_clump", .x)
}

#' @export
#' @rdname loci_ld_clump
loci_ld_clump.tbl_df <- function(.x, .col = "genotypes", ...) {
  stopifnot_gen_tibble(.x)
  # check n_cores are available
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_ld_clump only works with the genotypes column")
  }
  loci_ld_clump(.x$genotypes, ...)
}


#' @export
#' @rdname loci_ld_clump
loci_ld_clump.vctrs_bigSNP <- function(
    .x,
    .col = "genotypes",
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
    .old_opt <- getOption("bigstatsr.check.parallel.blas", TRUE)
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = .old_opt), add = TRUE)
  }

  # check that the loci have not been resorted
  # check that big_index in the loci table is an increasing sequence of indeces
  if (is.unsorted(show_loci(.x)$big_index, strictly = TRUE)) {
    stop(paste(
      "Your loci have been resorted; first save the new file backed",
      "matrix with `gt_update_backing_file()"
    ))
  }

  if (gt_has_imputed(.x) && gt_uses_imputed(.x) == FALSE) {
    # not uses_imputed
    gt_set_imputed(.x, set = TRUE)
    on.exit(gt_set_imputed(.x, set = FALSE), add = TRUE)
  }

  is_loci_table_ordered(.x, error_on_false = TRUE)

  # get the FBM
  geno_fbm <- attr(.x, "fbm") # nolint
  # rows (individuals) that we want to use
  if (use_positions) {
    .positions <- rep(NA, ncol(geno_fbm))
    .positions[show_loci(.x)$big_index] <- show_loci(.x)$position
  } else {
    .positions <- NULL
  }
  # create a chromosome vector (fill gaps between bigsnpr and show_loci)
  .chromosome <- rep(2147483647L, ncol(geno_fbm))
  .chromosome[show_loci(.x)$big_index] <-
    cast_chromosome_to_int(show_loci(.x)$chromosome)
  # now figure out if we have any snp which have already been removed
  # those will go into `exclude`
  loci_not_in_tibble <-
    seq_len(ncol(geno_fbm))[
      !seq_len(ncol(geno_fbm)) %in% .gt_fbm_cols(.x)
    ] # nolint
  exclude <- c(loci_not_in_tibble, .gt_fbm_cols(.x)[exclude])
  if (length(exclude) == 0) {
    exclude <- NULL
  }

  # Normalize S to FBM-wide vector if needed
  if (!is.null(S)) {
    fbm_n <- ncol(geno_fbm)
    visible_n <- length(.gt_fbm_cols(.x))
    if (length(S) == visible_n) { # nolint start
      S_full <- rep(NA_real_, fbm_n)
      S_full[.gt_fbm_cols(.x)] <- S
      S <- S_full
    } else if (length(S) != fbm_n) { # nolint end
      stop("Length of 'S' must equal length(.gt_fbm_cols(.x)) or ncol(FBM).")
    }
  }

  # as long as we have more than one individual
  snp_clump_ids <- bigsnpr::snp_clumping(
    G = geno_fbm,
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
