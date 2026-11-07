#' Compute basic population global statistics
#'
#' This function computes basic population global statistics, following the
#' notation in Nei 1987 (which in turn is based on Nei and Chesser 1983):
#' - observed heterozygosity ( \eqn{\hat{h}_o}, column header `Ho`)
#' - expected heterozygosity, also known as gene diversity
#' ( \eqn{\hat{h}_s}, `Hs`)
#' - total heterozygosity ( \eqn{\hat{h}_t}, `Ht`)
#' - genetic differentiation between subpopulations (\eqn{D_{st}}, `Dst`)
#' - corrected total population diversity (\eqn{h'_t}, `Htp`)
#' - corrected genetic differentiation between subpopulations
#' (\eqn{D'_{st}}, `Dstp`)
#' - \eqn{\hat{F}_{ST}} (column header, `Fst`)
#' - corrected \eqn{\hat{F'}_{ST}} (column header `Fstp`)
#' - \eqn{\hat{F}_{IS}} (column header, `Fis`)
#' - Jost's \eqn{\hat{D}} (column header, `Dest`)
#'
#' @details We use the notation of Nei 1987. That notation was for loci with
#'   \eqn{m} alleles, but in our case we only have two alleles, so `m=2`.
#' - Within population observed heterozygosity \eqn{\hat{h}_o} for a locus with
#'   \eqn{m} alleles is defined as:\cr
#'   \eqn{\hat{h}_o= 1-\sum_{k=1}^{s} \sum_{i=1}^{m}  \hat{X}_{kii}/s}\cr
#'   where\cr \eqn{\hat{X}_{kii}} represents the proportion of homozygote
#'   \eqn{i} in the sample for the \eqn{k}th population and\cr \eqn{s} the
#'   number of populations,\cr following equation 7.38 in Nei(1987) on
#'   pp.164.\cr
#'
#' - Within population expected heterozygosity (gene diversity)
#'   \eqn{\hat{h}_s} for a locus with \eqn{m} alleles is defined as:\cr
#'   \eqn{\hat{h}_s=(\tilde{n}/(\tilde{n}-1))[1-\sum_{i=1}^{m}\bar{\hat{x}_i^2}-\hat{h}_o/2\tilde{n}]}\cr #nolint
#'   where \cr \eqn{\tilde{n}=s/\sum_k 1/n_k} (i.e the harmonic mean of
#'   \eqn{n_k}) and\cr \eqn{\bar{\hat{x}_i^2}=\sum_k \hat{x}_{ki}^2/s}\cr
#'   following equation 7.39 in Nei(1987) on pp.164.
#'
#' - Total heterozygosity (total gene diversity)
#'   \eqn{\hat{h}_t} for a locus with \eqn{m} alleles is defined as:\cr
#'   \eqn{\hat{h}_t = 1-\sum_{i=1}^{m} \bar{\hat{x}_i^2} +
#'   \hat{h}_s/(\tilde{n}s) - \hat{h}_o/(2\tilde{n}s)}\cr where \cr
#'   \eqn{\hat{x}_i=\sum_k \hat{x}_{ki}/s}\cr following equation 7.40 in
#'   Nei(1987) on pp.164.\cr
#'
#' - The amount of gene diversity among samples \eqn{D_{ST}} is defined as:\cr
#'   \eqn{D_{ST} = \hat{h}_t - \hat{h}_s}\cr following the equation provided in
#'   the text at the top of page 165 in Nei(1987).
#'
#' - The corrected amount of gene diversity among samples
#'   \eqn{D'_{ST}} is defined as:\cr
#'   \eqn{D'_{ST} = (s/(s-1))D'_{ST}}\cr following the equation provided in the
#'   text at the top of page 165 in Nei(1987).
#'
#' - Total corrected heterozygosity (total gene diversity)
#'   \eqn{\hat{h}_t} is defined as:\cr
#'   \eqn{\hat{h'}_t = \hat{h}_s + D'_{ST}}\cr following the equation provided
#'   in the text at the top of page 165 in Nei(1987).
#' - \eqn{\hat{F}_{IS}} is defined as:\cr
#'   \eqn{\hat{F}_{IS} = 1 - \hat{h}_o/\hat{h}_s}\cr following equation 7.41 in
#'   Nei(1987) on pp.164.\cr
#'
#' - \eqn{\hat{F}_{ST}} is defined as:\cr
#'   \eqn{\hat{F}_{ST} = 1 - \hat{h}_s/\hat{h}_t = D_{ST}/\hat{h}_t}\cr
#'   following equation 7.43 in Nei(1987) on pp.165.\cr
#'
#' - \eqn{\hat{F'}_{ST}} is defined as:\cr
#'   \eqn{\hat{F'}_{ST} = D'_{ST}/\hat{h'}_t}\cr following the explanation
#'   provided in the text at the top of page 165 in Nei(1987).
#'
#' - Jost's \eqn{\hat{D}} is defined as:\cr
#'   \eqn{\hat{D} = (s/(s-1))((\hat{h'}_t-\hat{h}_s)/(1-\hat{h}_s))}\cr as
#'   defined by Jost(2008)
#'
#'   All these statistics are first computed by locus, and then averaged across
#'   loci (including any monomorphic locus) to obtain genome-wide values. The
#'   function uses the same algorithm as `hierfstat::basic.stats()` but is
#'   optimized for speed and memory usage.
#'
#' @references Nei M, Chesser R (1983) Estimation of fixation indexes and gene
#'   diversities. Annals of Human Genetics, 47, 253-259.
#'
#'   Nei M. (1987) Molecular
#'   Evolutionary Genetics. Columbia University Press, pp. 164-165.
#'
#'   Jost L
#'   (2008) GST and its relatives do not measure differentiation. Molecular
#'   Ecology, 17, 4015-4026.
#'
#' @param .x a [`gen_tibble`] (usually grouped, as obtained by using
#'   [dplyr::group_by()]; use on a single population will return a number of
#'   quantities as NA/NaN).
#'   If grouped by more than one variable, groups are labelled by the values
#'   of the grouping variables pasted with "." (e.g. `pop1.F`).
#' @param by_locus boolean, determining whether the statistics should be
#'   returned by locus(TRUE), or as a single genome wide value (FALSE, the
#'   default).
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @returns a tibble of population statistics, with populations as rows and
#'   statistics as columns
#' @export
#' @seealso [hierfstat::basic.stats()]
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # Compute population global statistics
#' example_gt %>% pop_global_stats()
#'
#' # To return by locus, set by_locus = TRUE
#' example_gt %>% pop_global_stats(by_locus = TRUE)
# this code is adapted from hierfstat::basic.stats by Jerome Goudet
pop_global_stats <- function(
  .x,
  by_locus = FALSE,
  n_cores = bigstatsr::nb_cores()
) {
  stopifnot_diploid(.x)
  # compute the statistics for each locus in a single pass over the
  # genotypes (see .pop_stats_dip_pseudo() for the formulae); we only keep
  # the statistics by locus if needed, otherwise we just get their means
  pop_stats <- .pop_stats_dip_pseudo(
    .x,
    by_locus_global = by_locus,
    n_cores = n_cores
  )
  stat_names <- c(
    "Ho",
    "Hs",
    "Ht",
    "Dst",
    "Htp",
    "Dstp",
    "Fst",
    "Fstp",
    "Fis",
    "Dest"
  )

  if (by_locus) {
    res <- data.frame(pop_stats$global)
    names(res) <- stat_names
    return(res)
  } else {
    # overall summaries: the means across loci of each statistic (ignoring
    # NA and infinite values)
    overall <- pop_stats$mean_global_finite
    overall[7] <- overall[4] / overall[3]
    overall[8] <- overall[6] / overall[5]
    overall[9] <- 1 - overall[1] / overall[2]
    overall[10] <- overall[6] / (1 - overall[2])
    names(overall) <- stat_names
    return(overall)
  }
}

#' Compute population statistics in a single pass over the genotypes
#'
#' Internal function used by [pop_global_stats()], [pop_het_exp()],
#' [pop_het_obs()] and [pop_fis()] for diploid (and pseudohaploid) data.
#' Loci are processed in blocks: the genotypes of each group are summarised,
#' and for each locus we compute, for each group, the observed
#' heterozygosity `Ho`, the gene diversity `Hs` (Nei 1987 eqn 7.39, see
#' [pop_het_exp()]) and `Fis = 1 - Ho/Hs`, as well as the 10 global
#' statistics returned by
#' [pop_global_stats()] (based on means across groups, see the formulae in
#' its documentation). Statistics by locus are only stored if requested,
#' otherwise only their means across loci are returned (so memory use does
#' not depend on the number of loci). Means across groups and across loci
#' are accumulated as in [rowMeans()] and [colMeans()] (with
#' `na.rm = TRUE`).
#'
#' @param .x a [`gen_tibble`], grouped or not (in which case all individuals
#'   are treated as a single population)
#' @param by_locus_ho,by_locus_hs,by_locus_fis booleans, should matrices of
#'   the group statistics by locus (loci as rows, groups as columns) be
#'   returned?
#' @param by_locus_global boolean, should a matrix of the global statistics
#'   by locus (loci as rows, 10 statistics as columns) be returned?
#' @param n_cores number of cores to be used
#' @returns a list with elements `ho`, `hs`, `fis` and `global` (matrices by
#'   locus, empty if not requested), `mean_ho`, `mean_hs` and `mean_fis`
#'   (means across loci for each group), `mean_global` (means across loci of
#'   the global statistics) and `mean_global_finite` (the same, also
#'   ignoring infinite values), and `group_names` (the names of the groups)
#' @keywords internal
#' @noRd
.pop_stats_dip_pseudo <- function(
  .x,
  by_locus_ho = FALSE,
  by_locus_hs = FALSE,
  by_locus_fis = FALSE,
  by_locus_global = FALSE,
  n_cores = bigstatsr::nb_cores()
) {
  # get the populations if it is a grouped gen_tibble
  if (inherits(.x, "grouped_df")) {
    groups <- .group_info(.x)
  } else {
    # create a dummy pop
    groups <- list(n = 1L, ids = rep(0L, nrow(.x)), labels = "pop")
  }

  pop_stats <- pop_stats_dip_pseudo_cpp(
    .gt_get_fbm(.x),
    rowInd = .gt_fbm_rows(.x),
    colInd = .gt_fbm_cols(.x),
    groupIds = groups$ids,
    ngroups = groups$n,
    ploidy = indiv_ploidy(.x),
    by_locus_ho = by_locus_ho,
    by_locus_hs = by_locus_hs,
    by_locus_fis = by_locus_fis,
    by_locus_global = by_locus_global,
    block_size = 1024,
    ncores = n_cores
  )
  pop_stats$group_names <- groups$labels
  pop_stats
}
