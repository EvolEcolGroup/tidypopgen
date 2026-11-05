#' Compute the joint site frequency spectrum of each pair of populations
#'
#' Compute the two-dimensional (joint) site frequency spectrum (SFS) of each
#' pair of populations, i.e. the number of loci with \eqn{i} copies of the
#' derived (or minor) allele in the first population and \eqn{j} in the
#' second. Missing genotypes are dealt with as in `easySFS`, by projecting
#' each population down to `n_proj` allele copies (see [pop_sfs()] for
#' details). Each spectrum only uses the loci that have at least `n_proj`
#' called allele copies in both its populations, so its marginal spectra are
#' not, in general, equal to the spectra computed by [pop_sfs()]. A folded
#' spectrum is folded on the total number of copies across the two
#' populations (as in `dadi`), so that it counts copies of the allele that
#' is the minor one across the pair.
#'
#' @inheritParams pop_sfs
#' @param .x a grouped [`gen_tibble`] (as obtained by using
#'   [dplyr::group_by()]), with at least two groups.
#' @param fastest the order of the cells of each spectrum in the tidy
#'   output: "last" (the default) keeps the order of the populations, so
#'   that the rows are sorted by the number of allele copies in the first
#'   population, then in the second (the number of copies in the second
#'   population changes fastest, as in the files written by `dadi` and
#'   `easySFS`); "first" makes the number of copies in the first population
#'   change fastest.
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns If `type = "tidy"`, a tibble with one row per cell of the joint
#'   spectrum of each pair of populations, with the columns:
#'   - the groups of the pair (e.g. `population_1` and `population_2`; one
#'   column per grouping variable for each population of the pair, if the
#'   data are grouped by more than one variable);
#'   - `allele_count_1` and `allele_count_2`: the number of copies of the
#'   derived (or minor) allele in the first and second population (see
#'   `allele_count` in [pop_sfs()]);
#'   - `n_sites`: the (projected) number of sites;
#'   - `masked`: whether the cell is masked (see [pop_sfs()]).
#'
#'   The projection of each population is stored in the attribute `n_proj`.
#'   The tibble has class `pairwise_sfs_tbl`, and can be plotted with
#'   [autoplot()][autoplot_sfs].
#'   If `type = "array"`, a list of `sfs` objects (one per pair, named with
#'   the labels of the two groups pasted with "."), each an array with the
#'   first population in the rows and the second in the columns (see
#'   [pop_sfs()]). Pairs are in the order given by [utils::combn()].
#' @seealso [pop_sfs()], [nwise_pop_sfs()], [pop_sfs_preview()]
#' @rdname pairwise_pop_sfs
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' example_gt %>% pairwise_pop_sfs(n_proj = 2)
#'
#' # as arrays
#' example_gt %>% pairwise_pop_sfs(n_proj = 2, type = "array")
pairwise_pop_sfs <- function(.x, ...) {
  UseMethod("pairwise_pop_sfs", .x)
}

#' @export
#' @rdname pairwise_pop_sfs
pairwise_pop_sfs.gen_tbl <- function(
  .x,
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  fastest = c("last", "first"),
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
    stop("at least 2 groups are needed to compute pairwise spectra")
  }
  n_proj <- .sfs_n_proj(n_proj, data,
    n_cores = n_cores,
    block_size = block_size
  )

  # all pairs of groups, in a single pass over the genotypes
  pairwise_combn <- utils::combn(data$n_groups, 2)
  spectra_def <- lapply(seq_len(ncol(pairwise_combn)), function(i) {
    pairwise_combn[, i]
  })
  configs <- .sfs_config_counts(data,
    spectra = spectra_def,
    min_n = n_proj,
    n_cores = n_cores,
    block_size = block_size
  )
  spectra <- .sfs_from_configs(configs, spectra_def, n_proj,
    folded = folded,
    n_cores = n_cores
  )
  names(spectra) <- .pair_labels(data$labels, pairwise_combn)
  if (type == "array") {
    return(spectra)
  }

  cells <- lapply(spectra, .sfs_cells,
    fastest = fastest,
    count_names = c("allele_count_1", "allele_count_2")
  )
  res <- .sfs_bind_pair_rows(cells, data$groups, pairwise_combn)
  .sfs_tidy_attributes(res, n_proj, folded, class = "pairwise_sfs_tbl")
}

#' Bind the cells of the spectra of pairs of groups into a tidy tibble
#'
#' @param cells a list of data frames, the cells of the spectrum of each pair
#' @param groups the group information, as returned by `.group_info()`
#' @param pairwise_combn a matrix with two rows, giving the indices of the
#'   two groups of each pair
#' @returns a tibble with the groups of the pair of each row, with one column
#'   per grouping variable for each of the two groups (e.g. population_1 and
#'   population_2), followed by the columns of `cells`
#' @keywords internal
#' @noRd
.sfs_bind_pair_rows <- function(cells, groups, pairwise_combn) {
  pair_idx <- rep(seq_along(cells), vapply(cells, nrow, integer(1)))
  group_1 <- groups$keys[pairwise_combn[1, pair_idx], , drop = FALSE]
  names(group_1) <- paste0(groups$vars, "_1")
  group_2 <- groups$keys[pairwise_combn[2, pair_idx], , drop = FALSE]
  names(group_2) <- paste0(groups$vars, "_2")
  cells <- dplyr::bind_rows(unname(cells))
  # e.g. a grouping variable called allele_count gives a column
  # allele_count_1
  .check_group_col_names(c(names(group_1), names(group_2)), names(cells))
  dplyr::bind_cols(group_1, group_2, cells)
}

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
pairwise_pop_sfs.data.frame <- function(.x, ...) {
  pairwise_pop_sfs(restore_genotypes_tbl(.x, "pairwise_pop_sfs"), ...)
}
