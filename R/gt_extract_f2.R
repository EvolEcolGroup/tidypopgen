#' Compute and store blocked f2 statistics for ADMIXTOOLS 2
#'
#' @references Maier R, Patterson N (2024). admixtools: Inferring demographic
#'   history from genetic data. R package version 2.0.4,
#'   https://github.com/uqrmaie1/admixtools.
#'
#' This function prepares data for various *ADMIXTOOLS 2* functions from the
#' package *ADMIXTOOLS 2*. It takes a [`gen_tibble`], computes allele
#' frequencies and blocked f2-statistics, and writes the results to `outdir`. It
#' is equivalent to `admixtools::extract_f2()`.
#' @param .x a [`gen_tibble`]
#' @param outdir Directory where data will be stored.
#' @param blgsize SNP block size in Morgan. Default is 0.05 (5 cM). If `blgsize`
#'   is 100 or greater, it will be interpreted as base pair distance rather than
#'   centimorgan distance.
#' @param maxmem Maximum amount of memory to be used. If the required amount of
#'   memory exceeds `maxmem`, allele frequency data will be split into blocks,
#'   and the computation will be performed separately on each block pair. This
#'   doesn't put a precise cap on the amount of memory used (it used to at some
#'   point). Set this parameter to lower values if you run out of memory while
#'   running this function. Set it to higher values if this function is too slow
#'   and you have lots of memory.
#' @param maxmiss Discard SNPs which are missing in a fraction of populations
#'   higher than `maxmiss`
#' @param minmaf Discard SNPs with minor allele frequency less than `minmaf`
#' @param maxmaf Discard SNPs with minor allele frequency greater than `maxmaf`
#' @param minac2 Discard SNPs with allele count lower than 2 in any population
#'   (default `FALSE`). This option should be set to `TRUE` when computing
#'   f3-statistics where one population consists mostly of pseudohaploid
#'   samples. Otherwise heterozygosity estimates and thus f3-estimates can be
#'   biased. `minac2 == 2` will discard SNPs with allele count lower than 2 in
#'   any non-singleton population (this option is experimental and is based on
#'   the hypothesis that using SNPs with allele count lower than 2 only leads to
#'   biases in non-singleton populations). Note that while the `minac2` option
#'   discards SNPs with allele count lower than 2 in any population, the
#'   \code{qp3pop} function will only discard SNPs with allele count lower than
#'   2 in the first (target) population (when the first argument is the prefix
#'   of a genotype file; i.e. it is applied directly to a genotype file, not via
#'   precomputing f2 from a [`gen_tibble`]).
#' @param outpop Keep only SNPs which are heterozygous in this population
#' @param outpop_scale Scale f2-statistics by the inverse `outpop`
#'   heterozygosity (`1/(p*(1-p))`). Providing `outpop` and setting
#'   `outpop_scale` to `TRUE` will give the same results as the original
#' *qpGraph* when the `outpop` parameter has been set, but it has the
#'   disadvantage of treating one population different from the others. This may
#'   limit the use of these f2-statistics for other models.
#' @param transitions Set this to `FALSE` to exclude transition SNPs
#' @param transversions Set this to `FALSE` to exclude transversion SNPs
#' @param overwrite Overwrite existing files in `outdir`
#' @param adjust_pseudohaploid Genotypes of pseudohaploid samples are usually
#'   coded as `0` or `2`, even though only one allele is observed.
#'   `adjust_pseudohaploid` ensures that the observed allele count increases
#'   only by `1` for each pseudohaploid sample. If `TRUE` (default), samples
#'   that don't have any genotypes coded as `1` among the first 1000 SNPs are
#'   automatically identified as pseudohaploid. This leads to slightly more
#'   accurate estimates of f-statistics. Setting this parameter to `FALSE`
#'   treats all samples as diploid and is equivalent to the *ADMIXTOOLS* `
#'   inbreed: NO` option. Setting `adjust_pseudohaploid` to an integer `n` will
#'   check the first `n` SNPs instead of the first 1000 SNPs. NOW DEPRECATED,
#'   set the ploidy of the `gen_tibble` with [gt_pseudohaploid()].
#' @param fst Write files with pairwise FST for every population pair. Setting
#'   this to FALSE can make `extract_f2` faster and will require less memory.
#' @param afprod  Write files with allele frequency products for every
#'   population pair. Setting this to FALSE can make `extract_f2` faster and
#'   will require less memory.
#' @param poly_only Specify whether SNPs with identical allele frequencies in
#'   every population should be discarded (`poly_only = TRUE`), or whether they
#'   should be used (`poly_only = FALSE`). By default (`poly_only = c("f2")`),
#'   these SNPs will be used to compute FST and allele frequency products, but
#'   not to compute f2 (this is the default option in the original ADMIXTOOLS).
#' @param apply_corr Apply small-sample-size correction when computing
#'   f2-statistics (default `TRUE`)
#' @param n_cores Parallelize computation across `n_cores` cores.
#' @param quiet Suppress printing of progress updates
#' @return SNP metadata (invisibly)
#' @export
#' @examplesIf rlang::is_installed("admixtools")
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' lobsters <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#' lobsters <- lobsters %>% group_by(population)
#' f2_path <- tempfile()
#' gt_extract_f2(lobsters, outdir = f2_path, quiet = TRUE)
#' admixtools::f2_from_precomp(f2_path, verbose = FALSE)
gt_extract_f2 <- function(
    .x,
    outdir = NULL,
    blgsize = 0.05,
    maxmem = 8000,
    maxmiss = 0,
    minmaf = 0,
    maxmaf = 0.5,
    minac2 = FALSE,
    outpop = NULL,
    outpop_scale = TRUE,
    transitions = TRUE,
    transversions = TRUE,
    overwrite = FALSE,
    adjust_pseudohaploid = NULL,
    fst = TRUE,
    afprod = TRUE,
    poly_only = c("f2"),
    apply_corr = TRUE,
    n_cores = 1,
    quiet = FALSE) {
  if (!requireNamespace("admixtools", quietly = TRUE)) {
    stop(
      "to use this function, first install package 'admixtools' with\n",
      "devtools::install_github('uqrmaie1/admixtools')"
    )
  }
  # deprecation of adjust_pseudohaploid
  if (!is.null(adjust_pseudohaploid)) {
    stop(
      "adjust_pseudohaploid is deprecated. Set the ploidy of the ",
      "`gen_tibble` with `gt_pseudohaploid()`"
    )
  }

  # parameters that don't make sense with a gen_tibble
  inds <- NULL
  pops <- NULL
  pops2 <- NULL
  auto_only <- FALSE
  keepsnps <- NULL

  if (!inherits(.x, "grouped_df")) {
    stop(".x should be a grouped df")
  }
  # if no outdir is given, create a subdirectory f2 in the path of the
  # gen_tibble rds
  if (is.null(outdir)) {
    outdir <- file.path(dirname(.gt_get_fbm(.x)$rds), "f2")
  }

  verbose <- !quiet
  afdat <- gt_to_aftable(
    .x,
    n_cores = n_cores
  )

  if (is.null(inds)) pops <- union(pops, pops2)

  afdat <- afdat %>% # nolint
    discard_from_aftable(
      maxmiss = maxmiss,
      minmaf = minmaf,
      maxmaf = maxmaf,
      minac2 = minac2,
      outpop = outpop,
      transitions = transitions,
      transversions = transversions,
      keepsnps = keepsnps,
      auto_only = auto_only,
      poly_only = FALSE
    )
  afdat$snpfile <- afdat$snpfile %>% mutate(poly = as.logical(cpp_is_polymorphic(afdat$afs))) # nolint
  if (sum(afdat$snpfile$poly) == 0) stop("There are no informative SNPs!")

  if (verbose) {
    message(paste0(
      nrow(afdat$afs),
      " SNPs remain after filtering. ",
      sum(afdat$snpfile$poly),
      " are polymorphic.\n"
    ))
  }

  if (isTRUE(poly_only)) poly_only <- c("f2", "ap", "fst")
  arrs <- afs_to_f2_blocks( # nolint
    afdat,
    outdir = outdir,
    overwrite = overwrite,
    maxmem = maxmem,
    poly_only = poly_only,
    pops1 = pops,
    pops2 = pops2,
    outpop = if (outpop_scale) outpop else NULL,
    blgsize = blgsize,
    afprod = afprod,
    fst = fst,
    apply_corr = apply_corr,
    n_cores = n_cores,
    verbose = verbose
  )

  if (verbose) message(paste0("Data written to ", outdir, "/\n"))
  invisible(afdat$snpfile)
}


#' Create an allele frequency table (aftable) for admixtools
#'
#' This function is equivalent to `anygeno_to_aftable()`, but the aftable is
#' created directly from the `gen_tibble`.
#' @param .x the [`gen_tibble`]
#' @param adjust_pseudohaploid Genotypes of pseudohaploid samples are
#' usually coded as `0` or `2`, even though only one allele is observed.
#' `adjust_pseudohaploid` ensures that the observed allele count increases
#' only by `1` for each pseudohaploid sample. If `TRUE` (default), samples
#' that don't have any genotypes coded as `1` among the first 1000 SNPs are
#' automatically identified as pseudohaploid. This leads to slightly more
#' accurate estimates of f-statistics. Setting this parameter to `FALSE`
#' treats all samples as diploid and is equivalent to the *ADMIXTOOLS* `
#' inbreed: NO` option. Setting `adjust_pseudohaploid` to an integer `n`
#' will check the first `n` SNPs instead of the first 1000 SNPs.
#' @param n_cores Parallelize computation across `n_cores` cores.
#' @returns a list of 3 elements: afs, counts, and snp
#' @keywords internal
#' @noRd
gt_to_aftable <- function(
    .x,
    adjust_pseudohaploid = TRUE,
    n_cores = bigstatsr::nb_cores()) {
  if (!inherits(.x, "grouped_df")) {
    stop(".x should be a grouped df")
  }
  geno_fbm <- .gt_get_fbm(.x)



  aftable <- grouped_alt_freq_dip_pseudo_cpp(
    BM = geno_fbm,
    rowInd = .gt_fbm_rows(.x),
    colInd = .gt_fbm_cols(.x),
    groupIds = dplyr::group_indices(.x) - 1,
    ngroups = max(dplyr::group_indices(.x)),
    ncores = n_cores,
    ploidy = indiv_ploidy(.x),
    as_counts = FALSE
  )

  aftable <- list(
    afs = aftable[, 1:((ncol(aftable) / 2))],
    counts = aftable[, ((ncol(aftable) / 2) + 1):(ncol(aftable))]
  )

  .group_levels <- .x %>%
    group_keys() %>%
    pull(1)
  dimnames(aftable$afs) <- list(loci_names(.x), .group_levels)
  dimnames(aftable$counts) <- list(loci_names(.x), .group_levels)

  loci_new_names <- c(
    SNP = "name",
    CHR = "chromosome",
    POS = "position",
    cm = "genetic_dist",
    A1 = "allele_alt",
    A2 = "allele_ref"
  )
  snp <- show_loci(.x) %>%
    select(-all_of("big_index")) %>%
    rename(dplyr::all_of(loci_new_names)) %>%
    dplyr::relocate(all_of(c("CHR", "SNP", "cm", "POS", "A1", "A2")))
  snp$A1[is.na(snp$A1)] <- "0"
  class(snp) <- c("spec_tbl_df", class(snp))
  # we are missing column specs, which are generated by readr::col_spec
  # do we need to make sure that CHR and SNP are character?
  snp$CHR <- as.character(snp$CHR)
  snp$SNP <- as.character(snp$SNP)

  aftable$snpfile <- snp
  return(aftable)
}


# import some functions that are not exported by admixtools
discard_from_aftable <- function(...) {
  utils::getFromNamespace("discard_from_aftable", "admixtools")(...)
}

afs_to_f2_blocks <- function(...) {
  utils::getFromNamespace("afs_to_f2_blocks", "admixtools")(...)
}

cpp_is_polymorphic <- function(...) {
  utils::getFromNamespace("cpp_is_polymorphic", "admixtools")(...)
}
