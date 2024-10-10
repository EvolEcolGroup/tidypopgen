#' Compute population specific FIS
#'
#' This function computes population specific FIS, using either the approach of Nei 1987 ((as computed by [hierfstat::basic.stats()]).) or of Weir and Goudet 2017
#' (as computed by [hierfstat::fis.dosage()]).
#' @references
#' Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press
#' Weir, BS and Goudet J (2017) A Unified Characterization of Population Structure and Relatedness. Genetics (2017) 206:2085
#' @param .x a grouped [`gen_tibble`] (as obtained by using [dplyr::group_by()])
#' @param method one of "Nei87" (based on Nei 1987) or "WG17" (for Weir and Goudet 2017) to compute FIS
#' @param by_locus boolean, determining whether FIS should be returned by locus(TRUE),
#' or as a single genome wide value (FALSE, the default). Note that this is only relevant for "Nei87",
#' as "WG17" always returns a single value.
#' @param include_global boolean determining whether, besides the population specific estiamtes, a global
#' estimate should be appended. Note that this will return a vector of n populations plus 1 (the global value),
#' or a matrix with n+1 columns if `by_locus=TRUE`.
#' @param allele_sharing_mat optional and only relevant for "WG17", the matrix of Allele Sharing returned by
#' [pairwise_allele_sharing()] with `as_matrix=TRUE`. As a number of statistics can be
#' derived from the Allele Sharing matrix,
#' it it sometimes more efficient to pre-compute this matrix.
#' @returns a vector of population specific fis (plus the global value if `include_global=TRUE`)
#' @export

pop_fis <- function(.x, method = c("Nei87", "WG17"), by_locus = FALSE, include_global=FALSE, allele_sharing_mat = NULL){
  method <- match.arg(method)
  if (method == "Nei87"){
    if (!is.null(allele_sharing_mat)){
      stop("allele_sharing_mat not relevant for Nei87")
    }
    pop_fis_nei87(.x, by_locus = by_locus, include_global = include_global)
  } else if (method == "WG17"){
    if (by_locus){
      stop("by_locus not implemented for WG17")
    }
    pop_fis_wg17(.x, include_global=include_global, allele_sharing_mat = allele_sharing_mat)
  }
}

pop_fis_nei87 <- function(.x, by_locus = FALSE, include_global=include_global,
                          n_cores = bigstatsr::nb_cores()){
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
    Fis = 1 - sHo/Hs

    colnames(Fis)<- .group_levels %>% dplyr::pull(1)
    if (by_locus){
      if (include_global){
        global <- (.x %>% pop_global_stats(by_locus = TRUE, n_cores = n_cores))$Fis
        Fis <- cbind(Fis, global)
      }
    } else {
      Fis <- colMeans(Fis, na.rm = TRUE)
      if (include_global){
         global <- (.x %>% pop_global_stats(by_locus = FALSE, n_cores = n_cores))["Fis"]
         names(global) <- "global"
         Fis <- c(Fis, global)
      }
    }
    return(Fis)
}


pop_fis_wg17 <- function(.x, include_global=FALSE, allele_sharing_mat = NULL){
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  if (is.null(allele_sharing_mat)){
    allele_sharing_mat <- pairwise_allele_sharing(.x, as_matrix = TRUE)
  }
  fis_by_pop <- hierfstat::fis.dosage(allele_sharing_mat,
             matching=TRUE,
             pop=group_indices(.x))
  names(fis_by_pop) <- c(dplyr::group_keys(.x) %>% pull(1),"global")
  if (include_global){
    return(fis_by_pop)
  } else {
    return(fis_by_pop[-length(fis_by_pop)])
  }
}
