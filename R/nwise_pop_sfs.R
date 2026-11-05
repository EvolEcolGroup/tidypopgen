#' Compute the joint site frequency spectrum of all populations
#'
#' Compute the multidimensional (joint) site frequency spectrum (SFS) of all
#' the populations of a grouped `gen_tibble`, i.e. the number of loci with
#' each combination of numbers of copies of the derived (or minor) allele
#' in the populations (the "multiSFS" of `easySFS` and `fastsimcoal2`).
#' Missing genotypes are dealt with as in `easySFS`, by projecting each
#' population down to `n_proj` allele copies (see [pop_sfs()] for details).
#' Only the loci with at least `n_proj` called allele copies in every
#' population are used. A folded spectrum is folded on the total number of
#' copies across all populations (as in `dadi`).
#'
#' Populations are in the order of the groups (see [dplyr::group_keys()]):
#' to change their order (which matters for `dadi` and `fastsimcoal2`),
#' change the order of the levels of the grouping variable (e.g. with
#' [factor()]); to use a subset of the populations, [filter()][dplyr::filter]
#' the `gen_tibble` (empty groups are dropped). The spectrum has
#' `prod(n_proj + 1)` cells, which grows quickly with the number of
#' populations.
#'
#' @inheritParams pairwise_pop_sfs
#' @param fastest the order of the cells in the tidy output: "last" (the
#'   default) keeps the order of the populations, so that the rows are sorted
#'   by the number of allele copies in the first population, then in the
#'   second, and so on (the number of copies in the last population changes
#'   fastest, as in the files written by `dadi` and `easySFS`); "first" makes
#'   the number of copies in the first population change fastest.
#' @param max_cells the maximum number of cells of the spectrum; an error is
#'   raised for larger spectra, to avoid running out of memory.
#' @returns If `type = "tidy"`, a tibble with one row per cell of the joint
#'   spectrum, with one column per population (named with the group labels)
#'   giving its number of copies of the derived (or minor) allele (see
#'   `allele_count` in [pop_sfs()]),
#'   `n_sites` (the projected number of sites) and `masked` (whether the cell
#'   is masked, see [pop_sfs()]). The projection of each population is
#'   stored in the attribute `n_proj`. If `type = "array"`, an `sfs` object
#'   (see [pop_sfs()]), an array with one dimension per population.
#' @seealso [pop_sfs()], [pairwise_pop_sfs()], [pop_sfs_preview()]
#' @rdname nwise_pop_sfs
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' example_gt %>% nwise_pop_sfs(n_proj = 2)
#'
#' # the first population changing fastest
#' example_gt %>% nwise_pop_sfs(n_proj = 2, fastest = "first")
nwise_pop_sfs <- function(.x, ...) {
  UseMethod("nwise_pop_sfs", .x)
}

#' @export
#' @rdname nwise_pop_sfs
nwise_pop_sfs.gen_tbl <- function(
  .x,
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  fastest = c("last", "first"),
  max_cells = 1e7,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
) {
  rlang::check_dots_empty()
  if (!inherits(.x, "gen_tbl") || !inherits(.x, "grouped_df")) {
    stop(".x should be a grouped gen_tibble")
  }
  type <- match.arg(type)
  fastest <- match.arg(fastest)
  .sfs_check_folded(folded)
  data <- .sfs_set_ancestral(.sfs_data(.x), ancestral, folded = folded)
  if (data$n_groups < 2) {
    stop("at least 2 groups are needed to compute a joint spectrum")
  }
  if (type == "tidy" && any(data$labels %in% c("n_sites", "masked"))) {
    stop("groups cannot be called 'n_sites' or 'masked' with type = 'tidy'")
  }
  n_proj <- .sfs_n_proj(n_proj, data,
    n_cores = n_cores,
    block_size = block_size
  )
  n_cells <- prod(as.numeric(n_proj) + 1)
  if (n_cells > max_cells) {
    stop(
      "the spectrum would have ", format(n_cells), " cells (more than ",
      "`max_cells`); use fewer populations or smaller projections, or ",
      "increase `max_cells`"
    )
  }
  configs <- .sfs_config_counts(data,
    spectra = list(seq_len(data$n_groups)),
    min_n = n_proj,
    n_cores = n_cores,
    block_size = block_size
  )
  fs <- .sfs_from_configs(configs, list(seq_len(data$n_groups)), n_proj,
    folded = folded,
    n_cores = n_cores
  )[[1]]
  if (type == "array") {
    return(fs)
  }
  res <- .sfs_cells(fs, fastest = fastest, count_names = data$labels)
  .sfs_tidy_attributes(res, n_proj, folded)
}

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
nwise_pop_sfs.data.frame <- function(.x, ...) {
  nwise_pop_sfs(restore_genotypes_tbl(.x, "nwise_pop_sfs"), ...)
}
