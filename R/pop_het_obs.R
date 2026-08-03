#' Compute the population observed heterozygosity
#'
#' This function computes population heterozygosity, using the formula of Nei
#' (1987).
#'
#' @details Within population observed heterozygosity \eqn{\hat{h}_o} for a
#'   locus with \eqn{m} alleles is defined as:\cr \eqn{\hat{h}_o=
#'   1-\sum_{k=1}^{s} \sum_{i=1}^{m}  \hat{X}_{kii}/s}\cr where\cr
#'   \eqn{\hat{X}_{kii}} represents the proportion of homozygote \eqn{i} in the
#'   sample for the \eqn{k}th population and\cr \eqn{s} the number of
#'   populations,\cr following equation 7.38 in Nei(1987) on pp.164. In our
#'   specific case, there are only two alleles, so \eqn{m=2}. For population
#'   specific estimates, the sum is done over a single value of \eqn{k}.
#'   \eqn{\hat{h}_o} at the genome level is simply the mean of the locus
#'   estimates.
#'
#'   This statistic is ploidy-agnostic: an individual is heterozygous at a
#'   locus whenever its dosage is neither zero nor its own ploidy, so
#'   `pop_het_obs()` works for diploid, polyploid, and populations that mix
#'   individuals of different ploidy. The genome-wide estimate when
#'   `include_global = TRUE` is the unweighted mean of the population
#'   estimates for polyploid/mixed-ploidy data; for diploid and pseudohaploid
#'   data it is computed via [pop_global_stats()], for consistency with
#'   previous releases.
#'
#' @references Nei M. (1987) Molecular Evolutionary Genetics. Columbia
#'   University Press
#'
#' @param .x a [`gen_tibble`] (usually grouped, as obtained by using
#'   [dplyr::group_by()], otherwise the full tibble will be considered as
#'   belonging to a single population).
#' @param by_locus boolean, determining whether Ho should be returned by
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
#' example_gt %>% pop_het_obs()
#'
#' # To include the global expected heterozygosity, set include_global = TRUE
#' example_gt %>% pop_het_obs(include_global = TRUE)
#'
#' # To return by locus, set by_locus = TRUE
#' example_gt %>% pop_het_obs(by_locus = TRUE)
# adapted from hierfstat
pop_het_obs <- function(
  .x,
  by_locus = FALSE,
  include_global = FALSE,
  n_cores = bigstatsr::nb_cores()
) {
  stopifnot_gen_tibble(.x)
  stopifnot_no_pseudohaploid(.x)
  # get the populations if it is a grouped gen_tibble
  if (inherits(.x, "grouped_df")) {
    .group_levels <- .x %>% group_keys()
    .group_ids <- dplyr::group_indices(.x) - 1
  } else {
    # create a dummy pop
    .group_levels <- tibble(population = "pop")
    .group_ids <- rep(0, nrow(.x))
  }

  # keep the diploid/pseudohaploid path on the original allele-frequency
  # kernel shared with pop_het_exp()/pop_fis()/pop_global_stats() (so
  # nothing changes, in behaviour or speed, for the common diploid case),
  # and only use the more general, ploidy-aware kernel when the data
  # actually contain polyploid or mixed-ploidy individuals.
  if (is_diploid_only(.x) || is_pseudohaploid(.x)) {
    pop_freqs_df <- grouped_summaries_dip_pseudo_cpp(
      .gt_get_fbm(.x),
      rowInd = .gt_fbm_rows(.x),
      colInd = .gt_fbm_cols(.x),
      groupIds = .group_ids,
      ngroups = nrow(.group_levels),
      ploidy = indiv_ploidy(.x),
      ncores = n_cores
    )
    Ho <- pop_freqs_df$het_obs # nolint
  } else {
    het_obs_df <- grouped_het_obs_polyploid(
      .gt_get_fbm(.x),
      rowInd = .gt_fbm_rows(.x),
      colInd = .gt_fbm_cols(.x),
      groupIds = .group_ids,
      ngroups = nrow(.group_levels),
      ploidy = indiv_ploidy(.x),
      ncores = n_cores
    )
    Ho <- het_obs_df$het_obs # nolint
  }
  colnames(Ho) <- .group_levels %>% dplyr::pull(1) # nolint
  if (include_global) {
    if (is_diploid_only(.x) || is_pseudohaploid(.x)) {
      global <-
        (.x %>% pop_global_stats(by_locus = TRUE, n_cores = n_cores))$Ho
    } else {
      # Hs, Fst, Fis etc. in pop_global_stats() rely on allele-frequency
      # formulae not yet generalised beyond diploid/pseudohaploid data, but
      # the global Ho (Nei 1987 eqn 7.38) is simply the unweighted mean of
      # the population estimates, which holds regardless of ploidy.
      global <- rowMeans(Ho, na.rm = TRUE)
    }
    Ho <- cbind(Ho, global) # nolint
  }

  if (by_locus) {
    return(Ho)
  } else {
    return(colMeans(Ho, na.rm = TRUE))
  }
}
