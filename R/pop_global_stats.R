#' Compute basic population global statistics
#'
#' This function computes basic population global statistics, following the notation in Nei 1987 (which in turn is based on Nei and Chesser 1983):
#' - observed heterozygosity (\eqn{h_o}, column header `Ho`)
#' - expected heterozygosity, also known as gene diversity (\eqn{h_s}, `Hs`)
#' - total heterozygosity (\eqn{h_t}, `Ht`)
#' - within population diversity (`Dst`)
#' - total population diversity (`Htp`)
#' - total population diversity (`Dstp`)
#' - Fst
#' - Fst prime
#' - Fis
#' - Dest
#' @details We use the notation of Nei 1987. That notation was for loci with \eqn{m} alleles, but in our case we only have two alleles, so `m=2`.
#' - Within population observed heterozygosity \eqn{\hat{h}_o} for a locus with \eqn{m} alleles is defined as:\cr
#'  \eqn{\hat{h}_o= 1-\sum_{k=1}^{s} \sum_{i=1}^{m}  \hat{X}_{kii}/s}\cr
#' where\cr
#' \eqn{\hat{X}_{kii}} represents the proportion of homozygote \eqn{i} in the sample for the \eqn{k}th population and\cr
#' \eqn{s} the number of populations,\cr
#' following equation 7.38 in Nei(1987) on pp.164.\cr
#'
#' - Within population expected heterozygosity (gene diversity)  \eqn{\hat{h}_s} for a locus with \eqn{m} alleles is defined as:\cr
#'  \eqn{\hat{h}_s=(\tilde{n}/(\tilde{n}-1))[1-\sum_{i=1}^{m}\bar{\hat{x}_i^2}-\hat{h}_o/2\tilde{n}]}\cr
#'  where \cr \eqn{\tilde{n}=s/\sum_k 1/n_k} (i.e the harmonic mean of \eqn{n_k}) and\cr
#'  \eqn{\bar{\hat{x}_i^2}=\sum_k \hat{x}_{ki}^2/s}\cr
#' following equation 7.39 in Nei(1987) on pp.164.
#'
#' - Total heterozygosity (total gene diversity) \eqn{\hat{h}_t} for a locus with \eqn{m} alleles is defined as:\cr
#' \eqn{\hat{h}_t = 1-\sum_{i=1}^{m} \bar{\hat{x}_i^2} + \hat{h}_s/(\tilde{n}s) - \hat{h}_o/(2\tilde{n}s)}\cr
#' where \cr
#'  \eqn{\hat{x}_i=\sum_k \hat{x}_{ki}/s}\cr
#' following equation 7.40 in Nei(1987) on pp.164.\cr
#' - The amount of gene diversity among samples \eqn{Dst} is defined as:\cr
#' \eqn{Dst=Ht - Hs}\cr
#'
#' - Total corrected heterozygosity (total gene diversity) \eqn{\hat{h}_t}
#'
#' @param .x a [`gen_tibble`] (usually grouped, as obtained by using [dplyr::group_by()]; use on
#' a single population will return a number of quantities as NA/NaN)
#' @param by_locus boolean, determining whether the statistics should be returned by locus(TRUE),
#' or as a single genome wide value (FALSE, the default).
#' @param n_cores number of cores to be used, it defaults to [bigstatsr::nb_cores()]
#' @returns a tibble of population statistics, with populations as rows and statistics as columns
#' @export


# this code is adapted from hierfstat::basic.stats by Jerome Goudet
pop_global_stats <- function(.x, by_locus = FALSE, n_cores = bigstatsr::nb_cores()){
  stopifnot_diploid(.x)
  # get the populations if it is a grouped gen_tibble
  if (inherits(.x,"grouped_df")){
    .group_levels <- .x %>% group_keys()
    .group_ids <- dplyr::group_indices(.x)-1
  } else { # create a dummy pop
    .group_levels = tibble(population="pop")
    .group_ids <- rep(0,nrow(.x))
  }

  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(.gt_get_bigsnp(.x)$genotypes,
                                       rowInd = .gt_bigsnp_rows(.x),
                                       colInd = .gt_bigsnp_cols(.x),
                                       groupIds = .group_ids,
                                       ngroups = nrow(.group_levels),
                                       ncores = n_cores)

  # get the number of individuals
  n <-pop_freqs_df$n/2
  # is this correct?
  sHo <-pop_freqs_df$het_obs
  mHo <- apply(sHo, 1, mean, na.rm = TRUE)

  # sum of squared frequencies
  sp2 <- pop_freqs_df$freq_alt^2+pop_freqs_df$freq_ref^2
#  Hs <- (1 - sp2 - sHo/2/n)
#  Hs <- n/(n - 1) * Hs
#  Fis = 1 - sHo/Hs

  np <- apply(n, 1, fun <- function(x) sum(!is.na(x)))
  # mean sample size over the populations
  mn <- apply(n, 1, fun <- function(x) {
    sum(!is.na(x))/sum(1/x[!is.na(x)])
  })
  # mean sum of square frequencies
  msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
  mp2 <- rowMeans(pop_freqs_df$freq_alt)^2+rowMeans(pop_freqs_df$freq_ref)^2
  mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
  Ht <- 1 - mp2 + mHs/mn/np - mHo/2/mn/np
  mFis = 1 - mHo/mHs

  Dst <- Ht - mHs
  Dstp <- np/(np - 1) * Dst
  Htp = mHs + Dstp
  Fst = Dst/Ht
  Fstp = Dstp/Htp
  Dest <- Dstp/(1 - mHs)
  res <- data.frame(cbind(mHo, mHs, Ht, Dst, Htp, Dstp, Fst,
                            Fstp, mFis, Dest))
  names(res) <- c("Ho", "Hs", "Ht", "Dst", "Htp", "Dstp",
                    "Fst", "Fstp", "Fis", "Dest")

  if (by_locus){
    return(res)
  } else {
    # overall summaries
    is.na(res) <- do.call(cbind, lapply(res, is.infinite))
    overall <- apply(res, 2, mean, na.rm = TRUE)
    overall[7] <- overall[4]/overall[3]
    overall[8] <- overall[6]/overall[5]
    overall[9] <- 1 - overall[1]/overall[2]
    overall[10] <- overall[6]/(1 - overall[2])
    names(overall) <- names(res)
    return(overall)
  }
}



