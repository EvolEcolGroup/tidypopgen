#' Compute population specific fis
#'
#' This function computes population specific fis (as computed by [hierfstat::fis.dosage()]).
#' @param .x a grouped [`gen_tibble`] (as obtained by using [dplyr::group_by()])
#' @param include_global boolean determing whether, besides the population specific fis, a global
#' fis should be appended. Note that this will return a vector of n populations plus 1 (the global value)
#' @param as_matrix optional, the matrix of Allele Sharing returned by
#' [pairwise_as()]. As a number of statistics can be derived from the AS matrix,
#' it it sometimes more efficient to precompute this matrix.
#' @returns a vector of population specific fis (plus the global value if `include_global=TRUE`)
#' @export

pop_fis <- function(.x, include_global=FALSE, as_matrix = NULL){
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  if (is.null(as_matrix)){
    as_matrix <- pairwise_as(.x)
  }
  fis_by_pop <- hierfstat::fis.dosage(as_matrix,
             matching=TRUE,
             pop=group_indices(.x))
  names(fis_by_pop) <- c(dplyr::group_keys(.x) %>% pull(1),"global")
  if (include_global){
    return(fis_by_pop)
  } else {
    return(fis_by_pop[-length(fis_by_pop)])
  }
}
