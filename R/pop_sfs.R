#' Compute the site frequency spectrum of each population
#'
#' Compute the one-dimensional site frequency spectrum (SFS) of each
#' population, i.e. the number of loci with 0, 1, 2, ... copies of the
#' alternate (or minor) allele. Missing genotypes are dealt with as in
#' `easySFS` (which relies on `dadi`), so that the spectra can be used for
#' demographic inference with `dadi`, `moments` or `fastsimcoal2`. See
#' [pairwise_pop_sfs()] for the joint spectra of pairs of populations, and
#' [nwise_pop_sfs()] for the joint spectrum of all populations.
#'
#' @details
#' **Missing data.** Because of missing genotypes, the number of allele
#' copies called at a locus varies from locus to locus. Each population is
#' therefore projected down to a fixed number of allele copies, `n_proj`
#' (Marth et al. 2004; Gutenkunst et al. 2009): a locus with \eqn{n} called
#' allele copies, \eqn{k} of which are alternate, contributes to cell
#' \eqn{j} of the spectrum the hypergeometric probability of drawing \eqn{j}
#' alternate copies when sampling `n_proj` of its \eqn{n} copies without
#' replacement. Loci with fewer than `n_proj` called copies are dropped. Each
#' locus is thus spread over several cells, and the spectrum contains
#' fractional numbers of sites. This is the approach of `easySFS` (with the
#' `-a` option, i.e. using all loci).
#'
#' **Choosing the projection.** Projecting to more allele copies gives more
#' resolution (more segregating sites, especially rare variants), but drops
#' more loci with missing data. Following the advice of the `dadi` manual and
#' `easySFS`, use [pop_sfs_preview()] to compute the number of segregating
#' sites for each possible projection, and choose a value that keeps most
#' segregating sites. If `n_proj` is not given, each population is projected
#' to the number of allele copies that maximises its number of segregating
#' sites (the `best` projection in [pop_sfs_preview()]).
#'
#' **Folding.** Unless `folded = FALSE`, the spectrum is folded (i.e. it
#' counts copies of the minor allele), as the reference allele is not
#' necessarily ancestral. If `folded = FALSE`, the reference allele is taken
#' to be ancestral, so that the spectrum counts copies of the alternate
#' (derived) allele: this is only meaningful if the alleles of the `gen_tibble`
#' are polarised (note that some operations, such as merging datasets with
#' [rbind()], can swap the reference and alternate alleles).
#'
#' **Masked cells.** As in `dadi`, the corners of the spectrum (sites
#' monomorphic in the projected sample) are masked, and so are the cells that
#' were folded out in a folded spectrum. They are reported in the `masked`
#' column (or the `mask` attribute), but keep their values.
#'
#' Each spectrum only uses the loci that have enough called allele copies in
#' its own populations. Diploid, polyploid, mixed-ploidy and pseudohaploid
#' data are supported: each individual contributes as many allele copies as
#' its ploidy (one for pseudohaploids). The spectrum is computed from the
#' genotypes currently in use: if the `gen_tibble` uses imputed genotypes (see
#' [gt_set_imputed()]), these are treated as called genotypes, so you will
#' usually want to use `gt_set_imputed(x, FALSE)` first.
#'
#' Note that `easySFS` samples a single SNP per `CHROM` (i.e. per RAD locus)
#' by default; to do the same, thin the loci (e.g. with [select_loci()])
#' before computing the spectra.
#'
#' @references Marth GT, Czabarka E, Murvai J, Sherry ST (2004) The allele
#'   frequency spectrum in genome-wide human variation data reveals signals of
#'   differential demographic history in three large world populations.
#'   Genetics 166: 351-372.
#'
#'   Gutenkunst RN, Hernandez RD, Williamson SH, Bustamante CD (2009)
#'   Inferring the joint demographic history of multiple populations from
#'   multidimensional SNP frequency data. PLoS Genetics 5: e1000695.
#'
#'   Overcast I. easySFS. \url{https://github.com/isaacovercast/easySFS}
#' @param .x a [`gen_tibble`] (grouped by population, or ungrouped to compute
#'   a single spectrum for all individuals), or a vector of class
#'   `vctrs_bigSNP` (usually the `genotypes` column of a [`gen_tibble`]).
#' @param n_proj the number of allele copies to project each population to:
#'   either a single value used for all populations, a vector with one value
#'   per population (in the order of the groups), or a vector named with the
#'   group labels. If `NULL` (the default), each population is projected to
#'   the number of allele copies that maximises its number of segregating
#'   sites (see [pop_sfs_preview()]).
#' @param folded boolean, whether to fold the spectrum (the default, as the
#'   ancestral allele is not usually known). If `FALSE`, the reference allele
#'   is taken to be the ancestral one.
#' @param type the type of object to return, either "tidy" (the default) or
#'   "array".
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param block_size maximum number of loci read at once.
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns If `type = "tidy"`, a tibble with one row per cell of the
#'   spectrum of each population, with the columns:
#'   - `group` (or one column per grouping variable, if the data are grouped
#'   by more than one variable), absent for ungrouped data;
#'   - `allele_count`: the number of copies of the alternate allele (for an
#'   unfolded spectrum) or of the minor allele (for a folded spectrum);
#'   - `n_sites`: the (projected) number of sites;
#'   - `masked`: whether the cell is masked (see details).
#'
#'   The projection of each population is stored in the attribute `n_proj`.
#'   If `type = "array"`, a list of `sfs` objects (one per population, named
#'   by the group labels), or a single `sfs` object for ungrouped data. An
#'   `sfs` object is an array of the numbers of sites, with attributes
#'   giving its `mask`, the population labels (`pop_ids`), `n_proj`, whether
#'   it is `folded`, and the number of loci used (`n_loci`) and dropped for
#'   lack of called allele copies (`n_loci_dropped`). Use [tidy()] to turn
#'   it into a tibble.
#' @rdname pop_sfs
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # Preview the number of segregating sites for each projection
#' example_gt %>% pop_sfs_preview()
#'
#' # Compute the folded spectra, projecting each population to 2 copies
#' example_gt %>% pop_sfs(n_proj = 2)
#'
#' # Project each population to a different number of copies
#' example_gt %>% pop_sfs(n_proj = c(pop1 = 4, pop2 = 2, pop3 = 2))
#'
#' # Unfolded spectra, as arrays
#' example_gt %>% pop_sfs(n_proj = 2, folded = FALSE, type = "array")
pop_sfs <- function(.x, n_proj = NULL, folded = TRUE,
                    type = c("tidy", "array"), n_cores, block_size, ...) {
  UseMethod("pop_sfs", .x)
}

#' @export
#' @rdname pop_sfs
pop_sfs.tbl_df <- function(
  .x,
  n_proj = NULL,
  folded = TRUE,
  type = c("tidy", "array"),
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
) {
  # TODO this is a hack to deal with the class being dropped when going
  # through group_map
  stopifnot_gen_tibble(.x)
  pop_sfs(.x$genotypes,
    n_proj = n_proj, folded = folded, type = type,
    n_cores = n_cores, block_size = block_size, ...
  )
}

#' @export
#' @rdname pop_sfs
pop_sfs.vctrs_bigSNP <- function(
  .x,
  n_proj = NULL,
  folded = TRUE,
  type = c("tidy", "array"),
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
) {
  rlang::check_dots_empty()
  .pop_sfs(.x,
    n_proj = n_proj, folded = folded, type = type,
    n_cores = n_cores, block_size = block_size
  )
}

#' @export
#' @rdname pop_sfs
pop_sfs.grouped_df <- function(
  .x,
  n_proj = NULL,
  folded = TRUE,
  type = c("tidy", "array"),
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
) {
  rlang::check_dots_empty()
  stopifnot_gen_tibble(.x)
  .pop_sfs(.x,
    n_proj = n_proj, folded = folded, type = type,
    n_cores = n_cores, block_size = block_size
  )
}

# the work of pop_sfs(), for a vctrs_bigSNP or a grouped gen_tibble
.pop_sfs <- function(.x, n_proj, folded, type, n_cores, block_size) {
  type <- match.arg(type, c("tidy", "array"))
  .sfs_check_folded(folded)
  data <- .sfs_data(.x)
  spectra_def <- as.list(seq_len(data$n_groups))
  if (is.null(n_proj)) {
    # count the configurations of all loci with at least one called copy,
    # and choose the projections from them: configurations with fewer
    # copies than the projection are then dropped when projecting, so we
    # only need a single pass over the genotypes
    configs <- .sfs_config_counts(data,
      spectra = spectra_def,
      min_n = rep(1L, data$n_groups),
      n_cores = n_cores,
      block_size = block_size
    )
    n_proj <- .sfs_n_proj(NULL, data,
      n_cores = n_cores,
      block_size = block_size,
      configs_1d = configs
    )
  } else {
    n_proj <- .sfs_check_n_proj(n_proj, data$max_n)
    configs <- .sfs_config_counts(data,
      spectra = spectra_def,
      min_n = n_proj,
      n_cores = n_cores,
      block_size = block_size
    )
  }
  spectra <- .sfs_from_configs(configs, spectra_def, n_proj,
    folded = folded,
    n_cores = n_cores
  )
  names(spectra) <- data$labels

  if (type == "array") {
    if (is.null(data$groups)) {
      return(spectra[[1]])
    }
    return(spectra)
  }
  cells <- lapply(spectra, .sfs_cells, count_names = "allele_count")
  if (is.null(data$groups)) {
    res <- cells[[1]]
  } else {
    res <- .bind_group_rows(cells, data$groups)
  }
  .sfs_tidy_attributes(res, n_proj, folded)
}

#' Preview the projections of the site frequency spectrum
#'
#' For each population, compute the number of segregating sites in its
#' site frequency spectrum projected to each possible number of allele
#' copies (from 2 to the number of allele copies of the population), as the
#' `--preview` option of `easySFS`. This is used to choose `n_proj` in
#' [pop_sfs()], [pairwise_pop_sfs()] and [nwise_pop_sfs()]: projecting to more
#' allele copies gives more resolution (more segregating sites, especially
#' rare variants), but drops more loci with missing data. The `dadi` manual
#' recommends maximising the number of segregating sites, while keeping a
#' reasonable sample size.
#'
#' The number of segregating sites for a projection to \eqn{p} copies is
#' the expected number of loci that are polymorphic in a sample of \eqn{p}
#' copies drawn without replacement from the called copies of each locus
#' (loci with fewer than \eqn{p} called copies are dropped), i.e. the sum of
#' the cells of the projected spectrum that are not corners. It is the same
#' for the folded and unfolded spectrum.
#'
#' @param .x a [`gen_tibble`], grouped by population (or ungrouped, to treat
#'   all individuals as a single population).
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param block_size maximum number of loci read at once.
#' @returns a tibble with columns `group` (or one column per grouping
#'   variable, if the data are grouped by more than one variable; absent for
#'   ungrouped data), `n_proj` (the number of allele copies of the
#'   projection), `n_seg_sites` (the number of segregating sites, not rounded)
#'   and `n_loci` (the number of loci with at least `n_proj` called allele
#'   copies), and `best`, which is `TRUE` for the projection that maximises
#'   the number of segregating sites (the largest one in case of ties), used
#'   by default by [pop_sfs()], [pairwise_pop_sfs()] and [nwise_pop_sfs()].
#' @seealso [pop_sfs()]
#' @export
#' @examples
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' example_gt %>% pop_sfs_preview()
#'
#' # the projection chosen for each population
#' example_gt %>%
#'   pop_sfs_preview() %>%
#'   dplyr::filter(best)
pop_sfs_preview <- function(
  .x,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1)
) {
  stopifnot_gen_tibble(.x)
  data <- .sfs_data(.x)
  configs <- .sfs_config_counts(data,
    spectra = as.list(seq_len(data$n_groups)),
    min_n = rep(1L, data$n_groups),
    n_cores = n_cores,
    block_size = block_size
  )
  previews <- .sfs_preview_from_configs(configs, data$max_n)
  if (is.null(data$groups)) {
    tibble::as_tibble(previews[[1]])
  } else {
    .bind_group_rows(previews, data$groups)
  }
}

# check the folded argument
.sfs_check_folded <- function(folded) {
  if (!is.logical(folded) || length(folded) != 1 || is.na(folded)) {
    stop("`folded` should be TRUE or FALSE")
  }
  invisible(TRUE)
}

# add the projections and folding to a tidy spectrum, as attributes
.sfs_tidy_attributes <- function(res, n_proj, folded) {
  attr(res, "n_proj") <- n_proj
  attr(res, "folded") <- folded
  res
}
