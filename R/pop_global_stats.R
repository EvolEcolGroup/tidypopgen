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
#'   quantities as NA/NaN)
#' @param by_locus boolean, determining whether the statistics should be
#'   returned by locus(TRUE), or as a single genome wide value (FALSE, the
#'   default).
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @returns a tibble of population statistics, with populations as rows and
#'   statistics as columns
#' @export
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

  # get the number of individuals
  n <- pop_freqs_df$n / 2
  # is this correct?
  sHo <- pop_freqs_df$het_obs # nolint
  mHo <- rowMeans(sHo, na.rm = TRUE) # nolint
  # sum of squared frequencies
  sp2 <- pop_freqs_df$freq_alt^2 + pop_freqs_df$freq_ref^2
  #  Hs <- (1 - sp2 - sHo/2/n) #nolint start
  #  Hs <- n/(n - 1) * Hs
  #  Fis = 1 - sHo/Hs #nolint end

  # np <- apply(n, 1, fun <- function(x) sum(!is.na(x))) # nolint start
  # # mean sample size over the populations
  # mn <- apply(
  #   n,
  #   1,
  #   fun <- function(x) {
  #     sum(!is.na(x)) / sum(1 / x[!is.na(x)])
  #   }
  # ) #nolint end

  np_mn <- compute_np_mn(n)
  np <- np_mn$np
  mn <- np_mn$mn
  # mean sum of square frequencies
  msp2 <- rowMeans(sp2, na.rm = TRUE) # nolint start
  mp2 <- rowMeans(pop_freqs_df$freq_alt)^2 + rowMeans(pop_freqs_df$freq_ref)^2
  mHs <- mn / (mn - 1) * (1 - msp2 - mHo / 2 / mn)
  Ht <- 1 - mp2 + mHs / mn / np - mHo / 2 / mn / np
  mFis <- 1 - mHo / mHs

  Dst <- Ht - mHs
  Dstp <- np / (np - 1) * Dst
  Htp <- mHs + Dstp
  Fst <- Dst / Ht
  Fstp <- Dstp / Htp
  Dest <- Dstp / (1 - mHs)
  res <- data.frame(cbind(
    mHo,
    mHs,
    Ht,
    Dst,
    Htp,
    Dstp,
    Fst,
    Fstp,
    mFis,
    Dest # nolint end
  ))
  names(res) <- c(
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
    return(res)
  } else {
    # overall summaries
    is.na(res) <- do.call(cbind, lapply(res, is.infinite))
    overall <- colMeans(res, na.rm = TRUE)
    overall[7] <- overall[4] / overall[3]
    overall[8] <- overall[6] / overall[5]
    overall[9] <- 1 - overall[1] / overall[2]
    overall[10] <- overall[6] / (1 - overall[2])
    names(overall) <- names(res)
    return(overall)
  }
}
