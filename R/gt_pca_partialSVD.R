#' PCA controlling for LD for `gen_tibble` objects
#'
#' This function performs Principal Component Analysis on a `gen_tibble`,
#' using a fast truncated SVD with initial pruning and then iterative removal
#' of long-range LD regions. This function is a wrapper for [bigsnpr::snp_autoSVD()]
#'
#' @param x a `gen_tbl` object
#' @param k Number of singular vectors/values to compute. Default is `10`.
#'   **This algorithm should be used to compute a few singular vectors/values.**
#' @param fun_scaling Usually this  can be left unset, as it defaults to
#' [bigsnpr::snp_scaleBinom()], which is the appropriate function for biallelic SNPs.
#' Alternatively it is possible to use  custom function
#' (see [bigsnpr::snp_autoSVD()] for details.
#'
#' @export


gt_pca_partialSVD <- function(x, k = 10, fun_scaling = bigsnpr::snp_scaleBinom()) {
  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  if (is.character(show_loci(x)$chromosome)){
    infos_chr <- as.numeric(factor(show_loci(x)$chromosome))
  } else {
    infos_chr <- show_loci(x)$chromosome
  }
  this_svd  <- bigstatsr::big_SVD(X$genotypes,
                                  k=k,
                                    ind.row = vctrs::vec_data(x$genotypes),
                                    ind.col = show_loci(x)$big_index,
                                    fun.scaling = fun_scaling,
                                    block.size= block_size(nrow(X$genotypes)))
  # add names to the scores (to match them to data later)
  rownames(this_svd$u)<-x$id
  this_svd$method <- "partialSVD"
  this_svd$call <- match.call()
  class(this_svd) <- c("gt_pca", class(this_svd))
  this_svd
}
