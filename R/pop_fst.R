#' Compute population specific Fst
#'
#' This function computes population specific Fst (as computed by [hierfstat::fst.dosage()]).
#' @param .x a grouped [`gen_tibble`] (as obtained by using [dplyr::group_by()])
#' @param include_global boolean determining whether, besides the population specific Fst, a global
#' Fst should be appended. Note that this will return a vector of n populations plus 1 (the global value)
#' @param allele_sharing_mat optional, the matrix of Allele Sharing returned by
#' [pairwise_allele_sharing()] with `as_matrix=TRUE`. As a number of statistics can be
#' derived from the Allele Sharing matrix,
#' @returns a vector of population specific Fst (plus the global value if `include_global=TRUE`)
#' @export

pop_fst <- function(.x, include_global=FALSE, allele_sharing_mat = NULL){
  if (!inherits(.x,"grouped_df")){
    stop (".x should be a grouped df")
  }
  if (is.null(allele_sharing_mat)){
    allele_sharing_mat <- pairwise_allele_sharing(.x, as_matrix = TRUE)
  }
  fst_by_pop <- hierfstat::fst.dosage(allele_sharing_mat,
             matching=TRUE,
             pop=group_indices(.x))
  names(fst_by_pop) <- c(dplyr::group_keys(.x) %>% pull(1),"global")
  if (include_global){
    return(fst_by_pop)
  } else {
    return(fst_by_pop[-length(fst_by_pop)])
  }
}
