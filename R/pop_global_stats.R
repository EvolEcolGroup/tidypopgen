#' Compute basic population global statistics
#'
#' This function computes basic population global statistics, following the formulas used in `hietfstat::basic.stats`.
#'
#' @details TO BE CORRECTED (THIS IS INCOMPLETE, AND FORMULAE WERE FILLED IN BY COPILOT)
#' The statistics computed are:
#' - observed heterozygosity (`Ho`)
#' - expected heterozygosity (`Hs`)
#' - total heterozygosity (`Ht`)
#' - within population diversity (`Dst`)
#' - total population diversity (`Htp`)
#' - total population diversity (`Dstp`)
#' - Fst
#' - Fst prime
#' - Fis
#' - Dest
#' @details The formulas used are:
#' \eqn{Ho= 1-\sum_k \sum_i  P_{k,ii}/n_p}
#' where \eqn{P_{k,ii}} represents the proportion of homozygote \eqn{i} in sample \eqn{k} and
#' \eqn{n_{p}} the number of samples.
#' \eqn{Hs=\tilde{n}/(\tilde{n}-1)[1-\sum_i\bar{p_i^2}-Ho/2\tilde{n}]}
#' where \eqn{\tilde{n}=np/\sum_k 1/n_k} and
#' \eqn{\bar{p_i^2}=\sum_k p_{ki}^2/np}
#' \eqn{Ht=1-\sum_k \bar{p_k^2} + \bar{Hs}/\tilde{n}/np - Ho/2/\tilde{n}/np}
#' \eqn{Dst=Ht - Hs}
#'
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



