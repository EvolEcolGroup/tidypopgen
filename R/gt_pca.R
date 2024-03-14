#' Principal Component Analysis for `gen_tibble` objects
#'
#' This function performs Principal Component Analysis on a `gen_tibble`.
#' It is currently just a wrapper around [bigsnpr::snp_autoSVD()]
#'
#' @param x a `gen_tbl` object
#' @param k the number of principal components to return/compute
#' @param method only "autoSVD" is currently implemented.
#' @param ... further parameters to pass to the appropriate method
#' @returns a `gt_pca` object, which is a subclass of `bigSVD`
#'
#' @export

## Look at to manipulate ellipses when passing arguments
#https://stackoverflow.com/questions/60338114/updating-values-of-three-dot-ellipsis-in-r

gt_pca <- function(x, k = 10, method = c("autoSVD","randomSVD"),...){
  method <- match.arg(method)
  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  # check that in the dots we don't have ind.col or ind.row
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  this_svd  <- bigsnpr::snp_autoSVD(X$genotypes,
                      k = k,
                      infos.chr = show_loci(x)$chromosome,
                      infos.pos = show_loci(x)$position,
                      ind.row = vctrs::vec_data(x$genotypes),
                      ind.col = show_loci(x)$big_index,
                      ...)
  # add names to the scores (to match them to data later)
  rownames(this_svd$u)<-x$id
  this_svd$method <- method
  this_svd$call <- match.call()
  class(this_svd) <- c("gt_pca",class(this_svd))
  this_svd
}


# a print method
#' @method print gt_pca
#' @export
print.gt_pca <- function(x, ...){
  cat(" === PCA of gen_tibble object ===")
  cat("\nMethod: ")
  print(x$method)
  cat("\nCall ($call):")
  print(x$call)
  cat("\nEigenvalues ($d):\n", round(utils::head(x$d,6),3), ifelse(length(x$d)>6, "...\n", "\n") )
  cat("\nPrincipal component scores ($u):\n matrix with", nrow(x$u), "rows (individuals) and", ncol(x$u), "columns (axes)", "\n")
  cat("\nLoadings (Principal axes) ($v):\n matrix with", nrow(x$v), "rows (SNPs) and", ncol(x$v), "columns (axes)", "\n")
  cat("\n")
}
