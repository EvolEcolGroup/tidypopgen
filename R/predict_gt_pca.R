#' Predict scores of a PCA
#'
#' Predict the PCA scores for a [`gt_pca`], either for the original data or
#' projecting new data.
#'
#' @param object the [`gt_pca`] object
#' @param new_data a gen_tibble if scores are requested for a new dataset
#' @param impute_to_center boolean on whether to impute missing values in
#' `new_data` to the mean values used
#' to center the pca. This option is used to e.g. project ancient data onto a PCA
#' fitted to modern data. It defaults to TRUE.
#' @param prediction_type a string taking the value of "simple" and/or OADP (Online Augmentation, Decomposition, and Procrustes (OADP) projection)
#' @param block_size number of loci read simultaneously (larger values will speed up
#' computation, but require more memory)
#' @param n_cores number of cores
#' @param ... no used
#' @returns a matrix of predictions, with samples as rows and components as columns. The number
#' of components depends on how many were estimated in the [`gt_pca`] object. If prediction
#' type is equal to c("simple","OADP"), then a list of two matrices is returned
#' @references Zhang et al (2020). Fast and robust ancestry prediction using
#' principal component analysis  36(11): 3439–3446.
#' @rdname predict_gt_pca
#' @export

# this is a modified version of bigstatsr::predict.big_SVD
predict.gt_pca <- function(object, new_data=NULL, impute_to_center = TRUE,
                           prediction_type = "simple",
                           block_size = NULL,
                           n_cores = 1, ...){
  rlang::check_dots_empty()
  if (!all(prediction_type %in% c("simple", "OADP"))){
    stop("prediction_type can only take values 'simple' or 'OADP'")
  }

  if (is.null(new_data)) {
    # U * D
    sweep(object$u, 2, object$d, '*')
  } else {
    if (!inherits(new_data,"gen_tbl")){
      stop ("new_data should be a gen_tibble")
    }
    # check the new_data have the same loci as the dataset used to build the pca
    if (!all(object$loci$name %in% show_loci(new_data)$name)){
      stop("loci used in object are not present in new_data")
    }
    # get id of loci in new_data
    loci_subset <- match(object$loci$name, show_loci(new_data)$name)
    if (!all(all(show_loci(new_data)$allele_ref[loci_subset]==object$loci$allele_ref),
             all(show_loci(new_data)$allele_alt[loci_subset]==object$loci$allele_alt))){
      stop("ref and alt alleles differ between new_data and the data used to create the pca object")
    }

    if (!impute_to_center){
      if (gt_has_imputed(new_data) && !gt_uses_imputed(new_data)){
        gt_set_imputed(new_data, set = TRUE)
        on.exit(gt_set_imputed(new_data, set = FALSE))
      }

      if (is.null(block_size)){
        block_size <- bigstatsr::block_size(nrow(new_data))
      }
      # X * V
      XV <- bigstatsr::big_prodMat(.gt_get_bigsnp(new_data)$genotypes,
                                   object$v,
                                   ind.row = .gt_bigsnp_rows(new_data),
                                   ind.col = .gt_bigsnp_cols(new_data)[loci_subset],
                                   block.size = block_size,
                                   center = object$center,
                                   scale  = object$scale)
      # if we use OADP, then we need to compute Xnorm
      if ("OADP" %in% prediction_type){
        stop ("OADP currently only implemented for when `impute_to_center = TRUE`")
      }
    } else {

      X_norm <- bigstatsr::FBM(nrow(new_data), 1, init = 0)
      XV     <- bigstatsr::FBM(nrow(new_data), ncol(object$u), init = 0)

      bigstatsr::big_parallelize(
        .gt_get_bigsnp(new_data)$genotypes,
        p.FUN = fbm256_part_prod,
        ind = seq_along(loci_subset),
        ncores = n_cores,
        ind.row = .gt_bigsnp_rows(new_data),
        ind.col = .gt_bigsnp_cols(new_data)[loci_subset], #info_snp$`_NUM_ID_`[keep],
        center = object$center,
        scale = object$scale,
        V = object$v,
        XV = XV,
        X_norm = X_norm
      )

      if ("OADP" %in% prediction_type){
        oadp_proj <- utils::getFromNamespace("OADP_proj", "bigsnpr")(XV, X_norm, object$d, ncores = n_cores)
      }
    }
    if (all(c("simple","OADP") %in% prediction_type)){
      return(    list(
        simple_proj = XV[, , drop = FALSE],
        OADP_proj   = oadp_proj
      ))
    } else if ("simple" %in% prediction_type) {
      return(XV[,, drop = FALSE])
    } else {
      return (oadp_proj)
    }
    return(XV[,, drop = FALSE])
  }
}

#######################################################################################
# a port of bigsnpr::part_prod to work on standard fb256 matrices

fbm256_part_prod <- function(X, ind, ind.row, ind.col, center, scale, V, XV, X_norm) {

  res <- fbm256_prod_and_rowSumsSq(
    BM = X,
    ind_row = ind.row,
    ind_col = ind.col[ind],
    center  = center[ind],
    scale   = scale[ind],
    V       = V[ind, , drop = FALSE]
  )

  bigstatsr::big_increment(XV,     res[[1]], use_lock = TRUE)
  bigstatsr::big_increment(X_norm, res[[2]], use_lock = TRUE)
}

