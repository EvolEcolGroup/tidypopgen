#' Checks if a `gen_tibble` has been imputed
#'
#' This function checks if a dataset has been imputed. Note that having
#' imputation does not mean that the imputed values are used.
#'
#' @param x a `gen_tibble`
#' @returns boolean TRUE or FALSE depending on whether the dataset has been
#'   imputed
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # The initial gen_tibble contains no imputed values
#' example_gt %>% gt_has_imputed()
#'
#' # Now impute the gen_tibble
#' example_gt <- example_gt %>% gt_impute_simple()
#'
#' # And we can check it has been imputed
#' example_gt %>% gt_has_imputed()
gt_has_imputed <- function(x) {
  if (inherits(x, "gen_tbl")) {
    x <- x$genotypes
  }
  !is.null(attr(x, "imputed", exact = TRUE))
}

#' Checks if a `gen_tibble` uses imputed data
#'
#' This function checks if a dataset uses imputed data. Note that it is possible
#' to have a dataset that has been imputed but it is currently not using
#' imputation.
#'
#' @param x a `gen_tibble`
#' @returns boolean TRUE or FALSE depending on whether the dataset is using the
#'   imputed values
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Impute the gen_tibble
#' example_gt <- example_gt %>% gt_impute_simple()
#'
#' # Check whether the gen_tibble uses imputed values
#' example_gt %>% gt_uses_imputed()
gt_uses_imputed <- function(x) {
  if (inherits(x, "gen_tbl")) {
    x <- x$genotypes
  }
  if (!gt_has_imputed(x)) {
    stop("this dataset does not have any imputed values to use!")
  }
  if (
    identical(attr(x, "fbm")$code256, bigsnpr::CODE_IMPUTE_PRED) || # nolint
      identical(attr(x, "fbm")$code256, bigsnpr::CODE_DOSAGE)
  ) {
    # nolint
    return(TRUE) # nolint
  } else {
    return(FALSE) # nolint
  }
}

#' Sets a `gen_tibble` to use imputed data
#'
#' This function sets or unsets the use of imputed data. For some analysis
#' that does not allow for missing data (e.g. PCA with
#' `gt_pca(method = "random")`, or [gt_pca_autoSVD()]), we have to use
#' imputation, but for other analysis it might be preferable to allow for
#' missing data (e.g. [gt_pca()], [loci_ld_clump()] and [loci_ld_prune()]
#' handle missing genotypes, and use imputed ones only if they are set).
#'
#' @param x a `gen_tibble`
#' @param set a boolean defining whether imputed data should be used
#' @returns the gen_tibble, invisibly
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Impute the gen_tibble
#' example_gt <- example_gt %>% gt_impute_simple()
#'
#' # Check whether the gen_tibble uses imputed values
#' example_gt %>% gt_uses_imputed()
#'
#' # Set the gen_tibble to use imputed values
#' example_gt %>% gt_set_imputed(TRUE)
#'
#' # And check that the gen_tibble uses imputed values again
#' example_gt %>% gt_uses_imputed()
gt_set_imputed <- function(x, set = NULL) {
  if (inherits(x, "gen_tbl")) {
    x <- x$genotypes
  }
  # phased genotypes use their own code256 (see ?phased_encoding), which
  # would be overwritten
  if (is_phased(x)) {
    stop("gt_set_imputed() does not work on phased genotypes", call. = FALSE)
  }
  if (!is.logical(set)) {
    stop(paste(
      "set should be either TRUE (to use imputed)",
      "or FALSE (to use the raw calls"
    ))
  }
  if (!gt_has_imputed(x)) {
    stop("this dataset does not have imputed values to use!")
  }
  if (set == FALSE) {
    attr(x, "fbm")$code256 <- bigsnpr::CODE_012
  } else {
    if (attr(x, "imputed", exact = TRUE) %in% c("simple", "xgboost")) {
      attr(x, "fbm")$code256 <- bigsnpr::CODE_IMPUTE_PRED
    }
  }
  return(invisible(x))
}
