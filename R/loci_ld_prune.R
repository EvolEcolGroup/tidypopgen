#' Prune loci based on a Linkage Disequilibrium threshold, as PLINK
#'
#' This function prunes loci in LD as `plink --indep-pairwise` (PLINK 1.9)
#' does, and returns the same loci. A window of `window_size` loci (or of
#' `window_size` kb, with `window_in_kb = TRUE`) moves along each chromosome
#' by `step_size` loci. Within each window, pairs of loci whose squared
#' correlation (r^2) is above `thr_r2` are found in order, and the locus with
#' the lower minor allele frequency of each pair is removed (the later one
#' when their frequencies are equal), until no such pair remains. Loci with
#' no variance (e.g. monomorphic loci) are removed, as in PLINK; but, also as
#' in PLINK, a locus alone on its chromosome, or alone in the first window of
#' its chromosome (with a window in kb), is kept even if it has no variance.
#' Clumping ([loci_ld_clump()]) gives a better spread of loci along the
#' chromosomes.
#'
#' Missing genotypes do not need to be imputed: as in PLINK, r^2 is the
#' squared correlation of the allele counts of two loci over the individuals
#' genotyped at both. Pairs of loci genotyped together in fewer than
#' `min_shared` individuals have no reliable r^2: by default, they are not
#' considered linked (as in PLINK, for which r^2 is then undefined); with
#' `insufficient = "remove"`, the locus with the lower minor allele
#' frequency is removed, as if they were linked. Imputed genotypes are used
#' only if the `gen_tibble` is set to use them (see [gt_set_imputed()]).
#'
#' Differences from PLINK: all individuals are used (PLINK uses only
#' founders); all chromosomes are treated as autosomes (PLINK treats the
#' sex chromosomes of males as haploid); and loci on chromosome "0"
#' (unplaced) are pruned as those of any other chromosome, whereas PLINK
#' ignores them and lists them in neither its *.prune.in* nor its
#' *.prune.out* file (so that `--extract` with *.prune.in* drops them, but
#' `--exclude` with *.prune.out* keeps them).
#'
#' @param .x a [`gen_tibble`] object
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#'   directly to the function). This defaults to "genotypes" and can only take
#'   that value. There is no need for the user to set it, but it is included to
#'   resolve certain tidyselect operations.
#' @param window_size the size of the window: a number of loci (at least 2),
#'   or, with `window_in_kb = TRUE`, a distance in kb.
#' @param step_size the number of loci by which the window moves (PLINK
#'   recommends 1 with windows in kb). With a window of a number of loci, it
#'   should not be larger than `window_size`.
#' @param thr_r2 the threshold of r^2 above which two loci are linked (at
#'   least 0, and less than 1).
#' @param window_in_kb boolean: is `window_size` a distance in kb (`TRUE`), or
#'   a number of loci (`FALSE`, the default)?
#' @param min_shared the minimum number of individuals genotyped at both loci
#'   for their r^2 to be used. Default is `2`, as in PLINK (with fewer, r^2 is
#'   undefined).
#' @param insufficient what to do with pairs of loci genotyped together in
#'   fewer than `min_shared` individuals: `"keep"` (the default) both loci,
#'   or `"remove"` the one with the lower minor allele frequency, as if they
#'   were linked.
#' @param n_cores number of cores to be used (to read the genotypes; the
#'   pruning itself runs on one core, as in PLINK)
#' @param return_id boolean on whether the id of SNPs to keep should be
#'   returned. It defaults to FALSE, which returns a vector of booleans (TRUE or
#'   FALSE)
#' @param ... currently not used.
#' @returns a boolean vector indicating whether each locus is kept (as the
#'   loci of `plink --indep-pairwise` in the `.prune.in` file), or a vector
#'   of the indices of the loci kept (with `return_id = TRUE`)
#' @references Chang C.C. et al. (2015) Second-generation PLINK: rising to the
#'   challenge of larger and richer datasets. GigaScience 4: 7.
#' @seealso [loci_ld_clump()] for clumping.
#' @rdname loci_ld_prune
#' @export
#' @examples
#' # (the example has missing genotypes, which do not need to be imputed)
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # as plink --indep-pairwise 4 1 0.2
#' example_gt %>% loci_ld_prune(window_size = 4, step_size = 1, thr_r2 = 0.2)
#' # the indices of the loci kept
#' example_gt %>%
#'   loci_ld_prune(window_size = 4, step_size = 1, thr_r2 = 0.2,
#'     return_id = TRUE
#'   )
loci_ld_prune <- function(.x, .col = "genotypes", ...) {
  UseMethod("loci_ld_prune", .x)
}

#' @export
#' @rdname loci_ld_prune
loci_ld_prune.tbl_df <- function(.x, .col = "genotypes", ...) {
  stopifnot_has_genotypes(.x)
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_ld_prune only works with the genotypes column")
  }
  loci_ld_prune(.x$genotypes, ...)
}

#' @export
#' @rdname loci_ld_prune
loci_ld_prune.vctrs_bigSNP <- function(
  .x,
  .col = "genotypes",
  window_size,
  step_size,
  thr_r2,
  window_in_kb = FALSE,
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
  # the parameters, checked as PLINK does
  is_count <- function(x) {
    is.numeric(x) && length(x) == 1 && !is.na(x) && x >= 1 && x == round(x)
  }
  if (!is_count(window_size) || (!window_in_kb && window_size < 2)) {
    stop(
      "'window_size' should be an integer of at least 2 (or at least 1 ",
      "with window_in_kb = TRUE)"
    )
  }
  if (!is_count(step_size)) {
    stop("'step_size' should be a positive integer")
  }
  # (PLINK 1.9 accepts a larger step, but its window buffers then overflow:
  # it allocates them for window_size loci, and adds step_size loci to an
  # emptied window, so its results are undefined)
  if (!window_in_kb && step_size > window_size) {
    stop("'step_size' should not be larger than 'window_size'")
  }
  if (window_in_kb && step_size > 1) {
    warning("'step_size' should be 1 when the window is in kb (as in PLINK)")
  }
  valid_thr <- is.numeric(thr_r2) && length(thr_r2) == 1 && !is.na(thr_r2) &&
    thr_r2 >= 0 && thr_r2 < 1
  if (!valid_thr) {
    stop("'thr_r2' should be at least 0, and less than 1")
  }

  # check that the loci have not been resorted
  if (is.unsorted(show_loci(.x)$big_index, strictly = TRUE)) {
    stop(paste(
      "Your loci have been resorted; first save the new file backed",
      "matrix with `gt_update_backing_file()"
    ))
  }
  is_loci_table_ordered(.x, error_on_false = TRUE)

  loci <- show_loci(.x)
  # the minor allele frequencies, of the genotyped individuals
  maf <- loci_maf(.x, n_cores = n_cores)
  # prune each chromosome (each in increasing position)
  by_chrom <- split(seq_len(nrow(loci)), loci$chromosome)
  keep_id <- lapply(by_chrom, function(idx) {
    keep <- ld_prune_cpp(
      BM = attr(.x, "fbm"),
      rowInd = vctrs::vec_data(.x),
      colInd = loci$big_index[idx],
      pos = as.double(loci$position[idx]),
      window_size = as.integer(window_size),
      window_kb = window_in_kb,
      step = as.integer(step_size),
      thr = thr_r2,
      maf = maf[idx],
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
    keep_bool <- rep(FALSE, nrow(loci))
    keep_bool[to_keep_id] <- TRUE
    keep_bool
  }
}
