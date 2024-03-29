#' Principal Component Analysis for `gen_tibble` objects
#'
#' There are a number of PCA methods available for `gen_tibble` objects. They
#' are mostly designed to work on very large datasets, so they only compute
#' a limited number of components. For smaller datasets, `gt_partialSVD` allows
#' the use of partial (truncated) SVD to fit the PCA; this method is suitable when
#' the number of individuals is much smaller than the number of loci. For larger
#' dataset, `gt_randomSVD` is more appropriate. Finally, there is a method
#' specifically designed for dealing with LD in large datasets, `gt_autoSVD`.
#' Whilst this is arguably the best option, it is somewhat data hungry, and so
#' only suitable for very large datasets (hundreds of individuals with several hundred
#' thousands markers, or larger).
#' @name gt_pca
NULL


#' Autoplots for `gt_pca` objects
#'
#' This function produces simple plots to quickly inspect an object. They are
#' not customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots. For `gt_pca`, the following types of plots are available:
#' - `screeplot`: a plot of the eigenvalues of the principal components (currently
#' it plots the singular value)
#' - `scores` a scatterplot of the scores of each individual on two principal
#' components (defined by `pc`)
#' - `loadings` a plot of loadings of all loci for a given component (chosen with `pc`)
#'
#' @param object an object of class `gt_pca`
#' @param type the type of plot (one of "screeplot", "scores" and "loadings")
#' @param k the principal components to be plotted: for scores, a pair of values
#' e.g. c(1,2); for `loadings` either one or more values.
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @rdname autoplot_gt_pca
#' @export
autoplot.gt_pca <- function(object,
                type=c("screeplot", "scores","loadings"),
                k = NULL, ...)
{
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (type== "screeplot") {
    plot(object, type="screeplot")
  } else if (type == "scores"){
    if (is.null(k)){
      k <- c(1,2)
    }
    if (length(k)!=2){
      stop("for 'scores' plots, 'pc' should be a pair of values, e.g. c(1,2)")
      }
    plot(object, type = "scores", scores = k)
  } else if (type == "loadings"){
    if (is.null(k)){
      k <- 1
    }
    plot(object, type = "loadings", loadings = k)
  }
}
