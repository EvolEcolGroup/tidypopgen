#' Clump loci based on a Linkage Disequilibrium threshold
#'
#' This function uses clumping to remove SNPs at high LD. When used with its
#' default options, clumping based on MAF is similar to standard pruning (as
#' done by PLINK with "--indep-pairwise (size+1) 1 thr.r2", but it results in a
#' better spread of SNPs over the chromosome; see [loci_ld_prune()] for
#' pruning as in PLINK. The algorithm is that of [bigsnpr::snp_clumping()]:
#' the loci of each chromosome are visited in order of importance (`S`), and a
#' locus is kept unless its squared correlation with a more important locus
#' that is kept, within the window, is above `thr_r2`. See
#' https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html for more
#' information on the differences between pruning and clumping.
#'
#' Missing genotypes do not need to be imputed: as in PLINK (`--r2`,
#' `--indep-pairwise`), the squared correlation (r^2) of the allele counts of
#' two loci is computed over the individuals genotyped at both. With complete
#' (or imputed) genotypes, the results are those of
#' [bigsnpr::snp_clumping()]. Pairs of loci genotyped together in fewer than
#' `min_shared` individuals have no reliable r^2: by default, both loci are
#' kept (as in PLINK, for which r^2 is then undefined); with `insufficient =
#' "remove"`, the less important locus is removed, as if they were linked.
#' The pairs of loci whose r^2 is undefined because a locus does not vary
#' over the individuals genotyped at both are not linked. Imputed genotypes
#' are used only if the `gen_tibble` is set to use them (see
#' [gt_set_imputed()]).
#'
#' r^2 is computed on genotypes packed as bits (64 individuals per word), as
#' in PLINK; each chromosome is packed at once, which takes 3 bits per
#' genotype.
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
#'   for thresholding with respect to p-values associated with `S`. The
#'   indices are those of the loci of `.x` (whole numbers between 1 and the
#'   number of loci); a logical vector with one value per locus (`TRUE` for
#'   the loci to exclude) also works. Other values are an error.
#' @param min_shared the minimum number of individuals genotyped at both loci
#'   for their r^2 to be used. Default is `2`, as in PLINK (with fewer, r^2 is
#'   undefined).
#' @param insufficient what to do with pairs of loci genotyped together in
#'   fewer than `min_shared` individuals: `"keep"` (the default) both loci,
#'   or `"remove"` the less important one, as if they were linked.
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
#' @seealso [bigsnpr::snp_clumping()], whose algorithm this function
#'   follows; [loci_ld_prune()] for pruning as in PLINK.
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' # (the example has missing genotypes, which do not need to be imputed)
#' example_gt <- load_example_gt("gen_tbl")
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
  stopifnot_has_genotypes(.x)
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
  min_shared = 2,
  insufficient = c("keep", "remove"),
  n_cores = 1,
  return_id = FALSE,
  ...
) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  insufficient <- match.arg(insufficient)
  stopifnot_min_shared(min_shared)
  # check the number of cores (the kernel does not use the BLAS, so two
  # levels of parallelism are not a concern)
  .old_opt <- options(bigstatsr.check.parallel.blas = FALSE)
  on.exit(options(.old_opt), add = TRUE)
  bigparallelr::assert_cores(n_cores)

  # check that the loci have not been resorted
  # check that big_index in the loci table is an increasing sequence of indeces
  if (is.unsorted(show_loci(.x)$big_index, strictly = TRUE)) {
    stop(paste(
      "Your loci have been resorted; first save the new file backed",
      "matrix with `gt_update_backing_file()"
    ))
  }
  is_loci_table_ordered(.x, error_on_false = TRUE)

  loci <- show_loci(.x)
  n_loci <- nrow(loci)
  # the loci to exclude: indices of loci (a logical vector with one value
  # per locus, which earlier versions accepted, is turned into indices);
  # other values (out of range, not whole numbers, names) would otherwise
  # exclude nothing, or other loci, without an error
  if (is.logical(exclude) && length(exclude) == n_loci && !anyNA(exclude)) {
    exclude <- which(exclude)
  }
  if (length(exclude) > 0) {
    valid_exclude <- is.numeric(exclude) && !anyNA(exclude) &&
      all(exclude >= 1 & exclude <= n_loci & exclude == round(exclude))
    if (!valid_exclude) {
      stop(
        "'exclude' should be a vector of indices of loci (whole numbers ",
        "between 1 and ", n_loci, "), or a logical vector with one value ",
        "per locus"
      )
    }
  }
  # the statistic to prioritise loci: by default, the minor allele
  # frequencies (of the genotyped individuals)
  if (is.null(S)) {
    S <- loci_maf(.x, n_cores = n_cores) # nolint
  } else if (length(S) == ncol(attr(.x, "fbm"))) {
    # (one value per column of the backing matrix)
    S <- S[loci$big_index] # nolint
  } else if (length(S) != n_loci) {
    stop("Length of 'S' must equal length(.gt_fbm_cols(.x)) or ncol(FBM).")
  }
  # the loci used, split by chromosome (each in increasing position)
  used <- rep(TRUE, n_loci)
  used[exclude] <- FALSE
  by_chrom <- split(which(used), loci$chromosome[used])

  # clump each chromosome
  keep_id <- lapply(by_chrom, function(idx) {
    if (use_positions) {
      pos <- as.double(loci$position[idx])
      size_pos <- size * 1000 # kb to bp
    } else {
      # (the window counts the loci used)
      pos <- as.double(seq_along(idx))
      size_pos <- size
    }
    # the loci in order of importance (as bigsnpr, ties in their order)
    ord <- order(S[idx], decreasing = TRUE)
    keep <- ld_clump_cpp(
      BM = attr(.x, "fbm"),
      rowInd = vctrs::vec_data(.x),
      colInd = loci$big_index[idx],
      ordInd = ord,
      rankInd = match(seq_along(ord), ord),
      pos = pos,
      size = size_pos,
      thr = thr_r2,
      min_shared = min_shared,
      remove_insufficient = insufficient == "remove",
      ncores = n_cores
    )
    idx[keep == 1]
  })
  to_keep_id <- sort(unlist(keep_id, use.names = FALSE))

  if (return_id) {
    to_keep_id
  } else {
    keep_bool <- rep(FALSE, n_loci)
    keep_bool[to_keep_id] <- TRUE
    keep_bool
  }
}

# Check `min_shared`, the minimum number of individuals genotyped at both loci
# for their r^2 to be used (in loci_ld_clump() and loci_ld_prune())
stopifnot_min_shared <- function(min_shared) {
  valid <- is.numeric(min_shared) && length(min_shared) == 1 &&
    !is.na(min_shared) && min_shared >= 2 && min_shared == round(min_shared)
  if (!valid) {
    stop("'min_shared' should be an integer of at least 2")
  }
}
