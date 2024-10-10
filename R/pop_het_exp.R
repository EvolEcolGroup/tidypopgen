#' Compute the population expected heterozygosity
#'
#' This function computes expected population heterozygosity (also
#' referred to as gene diversity, to avoid the potentially misleading use of the term "expected" in this context), using the formula of Nei (1987).
#'
#' @details Within population expected heterozygosity (gene diversity) is defined as:
#'  \eqn{Hs=\tilde{n}/(\tilde{n}-1)[1-\sum_i\bar{p_i^2}-Ho/2\tilde{n}],}
#'
#'  where \eqn{\tilde{n}=np/\sum_k 1/n_k} and
#'  \eqn{\bar{p_i^2}=\sum_k p_{ki}^2/np}
#' following Nei(1987) in pp.164--5
#'
#' @references Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press
#'
#' @param .x a [`gen_tibble`] (usually grouped, as obtained by using [dplyr::group_by()], otherwise the full tibble
#' will be considered as belonging to a single population).
#' @param by_locus boolean, determining whether Hs should be returned by locus(TRUE),
#' or as a single genome wide value (FALSE, the default).
#' @param include_global boolean determining whether, besides the population specific estiamtes, a global
#' estimate should be appended. Note that this will return a vector of n populations plus 1 (the global value),
#' or a matrix with n+1 columns if `by_locus=TRUE`.
#' @param n_cores number of cores to be used, it defaults to [bigstatsr::nb_cores()]
#' @returns a vector of mean population observed heterozygosities (if `by_locus=FALSE`), or a matrix of
#' estimates by locus (rows are loci, columns are populations, `by_locus=TRUE`)
#' @export


# adapted from hierfstat
pop_het_exp <- function(.x, by_locus = FALSE, include_global = FALSE, n_cores = bigstatsr::nb_cores()){
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
  sHo <-pop_freqs_df$het_obs
  mHo <- apply(sHo, 1, mean, na.rm = TRUE)
  n <-pop_freqs_df$n/2
  # sum of squared frequencies
  sp2 <- pop_freqs_df$freq_alt^2+pop_freqs_df$freq_ref^2
  Hs <- (1 - sp2 - sHo/2/n)
  Hs <- n/(n - 1) * Hs

  colnames(Hs)<- .group_levels %>% dplyr::pull(1)
  if (include_global){
    global <- (.x %>% pop_global_stats(by_locus = TRUE, n_cores = n_cores))$Hs
    Hs <- cbind(Hs, global)
  }

  if (by_locus){
    return(Hs)
  } else {
    return(colMeans(Hs, na.rm = TRUE))
  }
}


# alias for gene diversity
#' @export
#' @rdname pop_het_exp
pop_gene_div <- function (...){
  pop_het_exp(...)
}
