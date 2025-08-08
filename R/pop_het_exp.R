#' Compute the population expected heterozygosity
#'
#' This function computes expected population heterozygosity (also referred to
#' as gene diversity, to avoid the potentially misleading use of the term
#' "expected" in this context), using the formula of Nei (1987).
#'
#' @details Within population expected heterozygosity (gene diversity)
#'   \eqn{\hat{h}_s} for a locus with \eqn{m} alleles is defined as:\cr
#'   \eqn{\hat{h}_s=\tilde{n}/(\tilde{n}-1)[1-\sum_{i}^{m}\bar{\hat{x}_i^2}-\hat{h}_o/2\tilde{n}]}\cr #nolint
#'
#'   where \cr \eqn{\tilde{n}=s/\sum_k 1/n_k} (i.e the harmonic mean of
#'   \eqn{n_k}) and\cr \eqn{\bar{\hat{x}_i^2}=\sum_k \hat{x}_{ki}^2/s}\cr
#'   following equation 7.39 in Nei(1987) on pp.164. In our specific case, there
#'   are only two alleles, so \eqn{m=2}. \eqn{\hat{h}_s} at the genome level for
#'   each population is simply the mean of the locus estimates for each
#'   population.
#'
#' @references Nei M. (1987) Molecular Evolutionary Genetics. Columbia
#'   University Press
#'
#' @param .x a [`gen_tibble`] (usually grouped, as obtained by using
#'   [dplyr::group_by()], otherwise the full tibble will be considered as
#'   belonging to a single population).
#' @param by_locus boolean, determining whether Hs should be returned by
#'   locus(TRUE), or as a single genome wide value (FALSE, the default).
#' @param include_global boolean determining whether, besides the population
#'   specific estimates, a global estimate should be appended. Note that this
#'   will return a vector of n populations plus 1 (the global value), or a
#'   matrix with n+1 columns if `by_locus=TRUE`.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @returns a vector of mean population observed heterozygosities (if
#'   `by_locus=FALSE`), or a matrix of estimates by locus (rows are loci,
#'   columns are populations, `by_locus=TRUE`)
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # Compute expected heterozygosity
#' example_gt %>% pop_het_exp()
#'
#' # To include the global expected heterozygosity, set include_global = TRUE
#' example_gt %>% pop_het_exp(include_global = TRUE)
#'
#' # To return by locus, set by_locus = TRUE
#' example_gt %>% pop_het_exp(by_locus = TRUE)
# adapted from hierfstat
pop_het_exp <- function(
    .x,
    by_locus = FALSE,
    include_global = FALSE,
    n_cores = bigstatsr::nb_cores()) {
  stopifnot_diploid(.x)
  # get the populations if it is a grouped gen_tibble
  if (inherits(.x, "grouped_df")) {
    .group_levels <- .x %>% group_keys()
    .group_ids <- dplyr::group_indices(.x) - 1
  } else {
    # create a dummy pop
    .group_levels <- tibble(population = "pop")
    .group_ids <- rep(0, nrow(.x))
  }

  # summarise population frequencies
  pop_freqs_df <- grouped_summaries_dip_pseudo_cpp(
    .gt_get_bigsnp(.x)$genotypes,
    rowInd = .gt_bigsnp_rows(.x),
    colInd = .gt_bigsnp_cols(.x),
    groupIds = .group_ids,
    ngroups = nrow(.group_levels),
    ploidy = indiv_ploidy(.x),
    ncores = n_cores
  )
  sHo <- pop_freqs_df$het_obs # nolint start
  mHo <- rowMeans(sHo, na.rm = TRUE)
  n <- pop_freqs_df$n / 2
  # sum of squared frequencies
  sp2 <- pop_freqs_df$freq_alt^2 + pop_freqs_df$freq_ref^2
  Hs <- (1 - sp2 - sHo / 2 / n)
  Hs <- n / (n - 1) * Hs # nolint end

  colnames(Hs) <- .group_levels %>% dplyr::pull(1) # nolint
  if (include_global) {
    global <- (.x %>% pop_global_stats(by_locus = TRUE, n_cores = n_cores))$Hs
    Hs <- cbind(Hs, global) # nolint
  }

  if (by_locus) {
    return(Hs)
  } else {
    return(colMeans(Hs, na.rm = TRUE))
  }
}


# alias for gene diversity
#' @export
#' @rdname pop_het_exp
pop_gene_div <- function(
    .x,
    by_locus = FALSE,
    include_global = FALSE,
    n_cores = bigstatsr::nb_cores()) {
  pop_het_exp(
    .x = .x,
    by_locus = by_locus,
    include_global = include_global,
    n_cores = n_cores
  )
}
