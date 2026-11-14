#' Block bootstrap of site frequency spectra
#'
#' Compute bootstrap replicates of the site frequency spectra of
#' [pop_sfs()], [pairwise_pop_sfs()] and [nwise_pop_sfs()], by resampling
#' blocks of loci with replacement. Replicates are used to estimate the
#' uncertainty of statistics computed from the spectra, or of the parameters
#' of demographic models fitted with `dadi`, `moments` or `fastsimcoal2`
#' (see [write_sfs()] to write the replicates in their formats).
#'
#' @details
#' The loci are divided into non-overlapping blocks that never span more
#' than one chromosome, defined as in the `jackknife_*` functions (e.g.
#' [jackknife_pop_tajimas_d()]): blocks of `block_size` SNPs, base pairs or
#' centiMorgans (the latter requires the `genetic_dist` column of the loci
#' table). Blocks should be larger than the extent of linkage
#' disequilibrium, so that they are approximately independent. Each
#' replicate draws, with replacement, as many blocks as there are (as the
#' bootstrap of `dadi`, `Misc.bootstraps_from_dd_chunks()`). With
#' `block_size = 1` and `size_unit = "snp"`, single SNPs are resampled (as the
#' `--resample-vcf` option of `easySFS`), which ignores linkage.
#'
#' A spectrum is a sum over loci, so the spectrum of a replicate is the sum
#' of the spectra of the blocks it drew (counting each block as many times as
#' it was drawn). The genotypes are therefore read only once, to compute the
#' spectrum of each block. Missing data are dealt with by projection, as in
#' [pop_sfs()], with the same projections (`n_proj`) for all replicates: if
#' `n_proj` is `NULL`, they are chosen from all the loci with
#' [pop_sfs_preview()]. The blocks are drawn with the random number
#' generator of R, so use [set.seed()] to make the replicates reproducible.
#'
#' @inheritParams pop_sfs
#' @inheritParams pairwise_pop_sfs
#' @param .x a [`gen_tibble`], grouped by population (or, for
#'   `bootstrap_pop_sfs()`, ungrouped to compute the spectrum of all
#'   individuals).
#' @param n_bootstrap the number of bootstrap replicates.
#' @param block_size the size of the blocks of loci that are resampled.
#' @param size_unit the unit of `block_size`: "snp" (the default), "cM" or
#'   "bp".
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns If `type = "tidy"`, a tibble with the columns of the tidy output
#'   of [pop_sfs()], [pairwise_pop_sfs()] or [nwise_pop_sfs()], and a column
#'   `replicate` (after the group columns) with the number of the replicate.
#'   If `type = "array"`, a list with one element per replicate, each of
#'   which is what [pop_sfs()], [pairwise_pop_sfs()] or [nwise_pop_sfs()]
#'   return with `type = "array"` (lists of `sfs` objects, or an `sfs` object).
#'   In both cases, the attribute `bootstrap_blocks` is a matrix with the
#'   number of times each block (columns, numbered in the order in which they
#'   appear in the loci table) was drawn in each replicate (rows), and the
#'   tibble also has the attributes `n_proj` and `folded`.
#' @seealso [pop_sfs()], [write_sfs()], [jackknife_stats_generic()]
#' @rdname bootstrap_sfs
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' set.seed(1)
#' # resample blocks of 2 SNPs
#' boot <- example_gt %>%
#'   bootstrap_pop_sfs(n_bootstrap = 10, block_size = 2, n_proj = 2)
#' boot
#'
#' # the 95% interval of the number of sites in each cell
#' boot %>%
#'   dplyr::group_by(group, allele_count) %>%
#'   dplyr::summarise(
#'     lower = stats::quantile(n_sites, 0.025),
#'     upper = stats::quantile(n_sites, 0.975)
#'   )
#'
#' # the joint spectra of each pair, as arrays
#' boot_2d <- example_gt %>%
#'   bootstrap_pairwise_pop_sfs(
#'     n_bootstrap = 5, block_size = 2, n_proj = 2,
#'     type = "array"
#'   )
#' boot_2d[[1]]$pop1.pop2
bootstrap_pop_sfs <- function(.x, ...) {
  UseMethod("bootstrap_pop_sfs", .x)
}

#' @export
#' @rdname bootstrap_sfs
bootstrap_pop_sfs.gen_tbl <- function(
  .x,
  n_bootstrap = 100,
  block_size,
  size_unit = c("snp", "cM", "bp"),
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  n_cores = bigstatsr::nb_cores(),
  ...
) {
  rlang::check_dots_empty()
  stopifnot_gen_tibble(.x)
  type <- match.arg(type)
  size_unit <- match.arg(size_unit)
  .sfs_check_folded(folded)
  .sfs_check_bootstrap_args(n_bootstrap, block_size)
  data <- .sfs_set_ancestral(.sfs_data(.x), ancestral, folded = folded)
  n_proj <- .sfs_n_proj(n_proj, data,
    n_cores = n_cores,
    block_size = .sfs_read_block_size(data)
  )
  spectra_def <- as.list(seq_len(data$n_groups))
  boot <- .sfs_bootstrap(data, spectra_def, n_proj,
    folded = folded, n_bootstrap = n_bootstrap, block_size = block_size,
    size_unit = size_unit, n_cores = n_cores
  )
  projected <- boot$projected
  if (type == "array") {
    res <- lapply(seq_len(n_bootstrap), function(i_rep) {
      spectra <- lapply(projected, .sfs_replicate, i_rep = i_rep)
      if (is.null(data$groups)) {
        return(spectra[[1]])
      }
      stats::setNames(spectra, data$labels)
    })
    attr(res, "bootstrap_blocks") <- boot$weights
    return(res)
  }
  cells <- lapply(projected, .sfs_replicate_cells,
    count_names = "allele_count"
  )
  if (is.null(data$groups)) {
    res <- cells[[1]]
  } else {
    res <- .bind_group_rows(cells, data$groups)
  }
  res <- .sfs_tidy_attributes(res, n_proj, folded)
  attr(res, "bootstrap_blocks") <- boot$weights
  res
}

#' @rdname bootstrap_sfs
#' @export
bootstrap_pairwise_pop_sfs <- function(.x, ...) {
  UseMethod("bootstrap_pairwise_pop_sfs", .x)
}

#' @export
#' @rdname bootstrap_sfs
bootstrap_pairwise_pop_sfs.gen_tbl <- function(
  .x,
  n_bootstrap = 100,
  block_size,
  size_unit = c("snp", "cM", "bp"),
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  fastest = c("last", "first"),
  n_cores = bigstatsr::nb_cores(),
  ...
) {
  rlang::check_dots_empty()
  if (!inherits(.x, "gen_tbl") || !inherits(.x, "grouped_df")) {
    stop(".x should be a grouped gen_tibble")
  }
  type <- match.arg(type)
  size_unit <- match.arg(size_unit)
  fastest <- match.arg(fastest)
  .sfs_check_folded(folded)
  .sfs_check_bootstrap_args(n_bootstrap, block_size)
  data <- .sfs_set_ancestral(.sfs_data(.x), ancestral, folded = folded)
  if (data$n_groups < 2) {
    stop("at least 2 groups are needed to compute pairwise spectra")
  }
  n_proj <- .sfs_n_proj(n_proj, data,
    n_cores = n_cores,
    block_size = .sfs_read_block_size(data)
  )
  pairwise_combn <- utils::combn(data$n_groups, 2)
  spectra_def <- lapply(seq_len(ncol(pairwise_combn)), function(i) {
    pairwise_combn[, i]
  })
  boot <- .sfs_bootstrap(data, spectra_def, n_proj,
    folded = folded, n_bootstrap = n_bootstrap, block_size = block_size,
    size_unit = size_unit, n_cores = n_cores
  )
  projected <- boot$projected
  pair_names <- .pair_labels(data$labels, pairwise_combn)
  if (type == "array") {
    res <- lapply(seq_len(n_bootstrap), function(i_rep) {
      stats::setNames(
        lapply(projected, .sfs_replicate, i_rep = i_rep),
        pair_names
      )
    })
    attr(res, "bootstrap_blocks") <- boot$weights
    return(res)
  }
  cells <- lapply(projected, .sfs_replicate_cells,
    fastest = fastest,
    count_names = c("allele_count_1", "allele_count_2")
  )
  # the groups of the pair of each row (as in pairwise_pop_sfs())
  res <- .sfs_bind_pair_rows(cells, data$groups, pairwise_combn)
  res <- .sfs_tidy_attributes(res, n_proj, folded)
  attr(res, "bootstrap_blocks") <- boot$weights
  res
}

#' @rdname bootstrap_sfs
#' @param max_cells the maximum number of cells of all the replicates (the
#'   number of cells of the joint spectrum times `n_bootstrap`); an error is
#'   raised for more cells, to avoid running out of memory.
#' @export
bootstrap_nwise_pop_sfs <- function(.x, ...) {
  UseMethod("bootstrap_nwise_pop_sfs", .x)
}

#' @export
#' @rdname bootstrap_sfs
bootstrap_nwise_pop_sfs.gen_tbl <- function(
  .x,
  n_bootstrap = 100,
  block_size,
  size_unit = c("snp", "cM", "bp"),
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  fastest = c("last", "first"),
  max_cells = 1e7,
  n_cores = bigstatsr::nb_cores(),
  ...
) {
  rlang::check_dots_empty()
  if (!inherits(.x, "gen_tbl") || !inherits(.x, "grouped_df")) {
    stop(".x should be a grouped gen_tibble")
  }
  type <- match.arg(type)
  size_unit <- match.arg(size_unit)
  fastest <- match.arg(fastest)
  .sfs_check_folded(folded)
  .sfs_check_bootstrap_args(n_bootstrap, block_size)
  data <- .sfs_set_ancestral(.sfs_data(.x), ancestral, folded = folded)
  if (data$n_groups < 2) {
    stop("at least 2 groups are needed to compute a joint spectrum")
  }
  reserved <- c("replicate", "n_sites", "masked")
  if (type == "tidy" && any(data$labels %in% reserved)) {
    stop(
      "groups cannot be called 'replicate', 'n_sites' or 'masked' with ",
      "type = 'tidy'"
    )
  }
  n_proj <- .sfs_n_proj(n_proj, data,
    n_cores = n_cores,
    block_size = .sfs_read_block_size(data)
  )
  n_cells <- prod(as.numeric(n_proj) + 1) * n_bootstrap
  if (n_cells > max_cells) {
    stop(
      "the replicates would have ", format(n_cells), " cells (more than ",
      "`max_cells`); use fewer replicates, fewer populations or smaller ",
      "projections, or increase `max_cells`"
    )
  }
  boot <- .sfs_bootstrap(data, list(seq_len(data$n_groups)), n_proj,
    folded = folded, n_bootstrap = n_bootstrap, block_size = block_size,
    size_unit = size_unit, n_cores = n_cores
  )
  projected <- boot$projected[[1]]
  if (type == "array") {
    res <- lapply(seq_len(n_bootstrap), function(i_rep) {
      .sfs_replicate(projected, i_rep = i_rep)
    })
    attr(res, "bootstrap_blocks") <- boot$weights
    return(res)
  }
  res <- .sfs_replicate_cells(projected,
    fastest = fastest,
    count_names = data$labels
  )
  res <- .sfs_tidy_attributes(res, n_proj, folded)
  attr(res, "bootstrap_blocks") <- boot$weights
  res
}

# check the number of bootstrap replicates and the size of the blocks
.sfs_check_bootstrap_args <- function(n_bootstrap, block_size) {
  is_count <- is.numeric(n_bootstrap) && length(n_bootstrap) == 1 &&
    !is.na(n_bootstrap) && n_bootstrap >= 1 &&
    n_bootstrap == round(n_bootstrap)
  if (!is_count) {
    stop("`n_bootstrap` should be a positive whole number")
  }
  if (missing(block_size) || is.null(block_size)) {
    stop("`block_size` is needed for the block bootstrap")
  }
  is_size <- is.numeric(block_size) && length(block_size) == 1 &&
    !is.na(block_size) && block_size > 0
  if (!is_size) {
    stop("`block_size` should be a positive number")
  }
  invisible(TRUE)
}

# the number of loci read at once by the kernels (as the default of the
# functions that compute the spectra)
.sfs_read_block_size <- function(data) {
  bigstatsr::block_size(length(data$genotypes), 1)
}

#' Compute bootstrap replicates of a set of spectra
#'
#' The configurations of each spectrum are counted in each block, in a single
#' pass over the genotypes, and each replicate is projected as the weighted
#' sum of its blocks.
#'
#' @param data the output of `.sfs_data()` (and `.sfs_set_ancestral()`)
#' @param spectra a list of integer vectors, each with the (1-based) groups of
#'   a spectrum
#' @param n_proj the projection of each group
#' @returns a list with `projected` (as returned by `.sfs_project()`, with one
#'   column per replicate) and `weights` (the number of times each block was
#'   drawn in each replicate)
#' @keywords internal
#' @noRd
.sfs_bootstrap <- function(data, spectra, n_proj, folded, n_bootstrap,
                           block_size, size_unit, n_cores) {
  blocks <- bootstrap_blocks(data$loci,
    n_bootstrap = n_bootstrap,
    block_size = block_size, size_unit = size_unit
  )
  configs <- .sfs_config_counts(data,
    spectra = spectra,
    min_n = n_proj,
    n_cores = n_cores,
    block_size = .sfs_read_block_size(data),
    block_ids = blocks$block_id
  )
  projected <- .sfs_project(configs, spectra, n_proj,
    folded = folded,
    n_cores = n_cores,
    weights = blocks$weights
  )
  list(projected = projected, weights = blocks$weights)
}

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
bootstrap_pop_sfs.data.frame <- function(.x, ...) {
  bootstrap_pop_sfs(restore_genotypes_tbl(.x, "bootstrap_pop_sfs"), ...)
}

#' @export
#' @noRd
bootstrap_pairwise_pop_sfs.data.frame <- function(.x, ...) {
  bootstrap_pairwise_pop_sfs(
    restore_genotypes_tbl(.x, "bootstrap_pairwise_pop_sfs"),
    ...
  )
}

#' @export
#' @noRd
bootstrap_nwise_pop_sfs.data.frame <- function(.x, ...) {
  bootstrap_nwise_pop_sfs(
    restore_genotypes_tbl(.x, "bootstrap_nwise_pop_sfs"),
    ...
  )
}
