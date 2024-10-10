#' Compute the population observed heterozygosity
#'
#' This function computes population heterozygosity, using the formula of Nei (1987).
#'
#' @details Within population observed heterozygosity is defined as:
#'  \eqn{Ho= 1-\sum_k \sum_i  P_{k,ii}/n_p}
#' where \eqn{P_{k,ii}} represents the proportion of homozygote \eqn{i} in sample \eqn{k} and
#' \eqn{n_{p}} the number of samples.
#' following Nei(1987) in pp.164--5
#'
#' @references Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press
#'
#' @param .x a [`gen_tibble`] (usually grouped, as obtained by using [dplyr::group_by()], otherwise the full tibble
#' will be considered as belonging to a single population).
#' @param by_locus boolean, determining whether Ho should be returned by locus(TRUE),
#' or as a single genome wide value (FALSE, the default).
#' @param include_global boolean determining whether, besides the population specific estiamtes, a global
#' estimate should be appended. Note that this will return a vector of n populations plus 1 (the global value),
#' or a matrix with n+1 columns if `by_locus=TRUE`.
#' @param n_cores number of cores to be used, it defaults to [bigstatsr::nb_cores()]
#' @returns a vector of mean population observed heterozygosities (if `by_locus=FALSE`), or a matrix of
#' estimates by locus (rows are loci, clumns are populations, `by_locus=TRUE`)
#' @export


# adapted from hierfstat
pop_het_obs <- function(.x, by_locus = FALSE, include_global = FALSE, n_cores = bigstatsr::nb_cores()){
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

  Ho <- pop_freqs_df$het_obs
  colnames(Ho)<- .group_levels %>% dplyr::pull(1)
  if (include_global){
    global <- (.x %>% pop_global_stats(by_locus = TRUE, n_cores = n_cores))$Ho
    Ho <- cbind(Ho, global)
  }

  if (by_locus){
    return(Ho)
  } else {
    return(colMeans(Ho, na.rm = TRUE))
  }
}
