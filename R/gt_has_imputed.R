#' Checks if a `gen_tibble` has been imputed
#'
#' This function checks if a dataset has been imputed. Note that having imputation
#' does not mean that the imputed values are used.
#'
#' @param x a `gen_tibble`
#' @returns boolean TRUE or FALSE depending on whether the dataset has been imputed
#' @export

gt_has_imputed <- function (x){
  !is.null(attr(x$genotypes,"imputed"))
}

#' Checks if a `gen_tibble` uses imputed data
#'
#' This function checks if a dataset uses imputed data. Note that it is possible
#' to have a dataset that has been imputed but it is currently not using imputation.
#'
#' @param x a `gen_tibble`
#' @returns boolean TRUE or FALSE depending on whether the dataset is using the imputed values
#' @export

gt_uses_imputed <- function (x){
  if (!gt_has_imputed(x)){
    stop("this dataset does not have any imputated values to use!")
  }
  if (identical(attr(x$genotypes,"bigsnp")$genotypes$code256, bigsnpr::CODE_IMPUTE_PRED) |
      identical(attr(x$genotypes,"bigsnp")$genotypes$code256, bigsnpr::CODE_DOSAGE)){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

#' Sets a `gen_tibble` to use imputed data
#'
#' This function sets or unsets the use of imputed data. For some analysis,
#' such as PCA, that does not allow for missing data, we have to use imputation,
#' but for other analysis it might be preferabble to allow for missing data.
#'
#' @param x a `gen_tibble`
#' @param set a boolean defining whether imputed data should be used
#' @returns the `gen_tibble` with the code switched to either use or
#' not use imputed values (depending on `set`)
#' @export

gt_set_imputed <- function (x, set = NULL){
  if (is.null(set)){
    stop ("set should be either TRUE (to use imputed) or FALSE (to use the raw calls")
  }
  if (!gt_has_imputed(x)){
    stop("this dataset does not have imputed values to use!")
  }
  if (set==FALSE){
    attr(x$genotypes,"bigsnp")$genotypes$code256 <- bigsnpr::CODE_012
  } else {
    if (attr(x$genotypes,"imputed") %in% c("simple","xgboost")){
      attr(x$genotypes,"bigsnp")$genotypes$code256 <- bigsnpr::CODE_IMPUTE_PRED
    }
  }
  x
}