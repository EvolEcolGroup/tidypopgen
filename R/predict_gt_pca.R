#' Predict scores of a PCA
#'
#' Predict the PCA scores for a [`gt_pca`], either for the original data or for new data.
#'
#' @param object the [`gt_pca`] object
#' @param new_data a gen_tibble if scores are requested for a new dataset
#' @param block_size number of loci read simultaneously (larger values will speed up
#' computation, but require more memory)
#' @param ... no used
#' @returns a matrix of predictions, with samples as rows and components as columns. The number
#' of components depends on how many were estimated in the [`gt_pca`] object
#' @rdname predict_gt_pca
#' @export

# this is a modified version of bigstatsr::predict.big_SVD
predict.gt_pca <- function(object, new_data=NULL,block_size = NULL, ...){
  rlang::check_dots_empty()
  if (is.null(new_data)) {
    # U * D
    sweep(object$u, 2, object$d, '*')
  } else {
    if (!inherits(new_data,"gen_tbl")){
      stop ("new_data should be a gen_tibble")
    }
    if (gt_has_imputed(new_data) && !gt_uses_imputed(new_data)){
      gt_set_imputed(new_data, set = TRUE)
      on.exit(gt_set_imputed(new_data, set = FALSE))
    }
    if (is.null(block_size)){
      block_size <- bigstatsr::block_size(nrow(new_data))
    }
    # X * V
    bigstatsr::big_prodMat(.gt_get_bigsnp(new_data)$genotypes,
                           object$v,
                ind.row = .gt_bigsnp_rows(new_data),
                ind.col = .gt_bigsnp_cols(new_data),
                block.size = block_size,
                center = object$center,
                scale  = object$scale)
  }
}