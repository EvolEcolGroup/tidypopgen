#' Autoplots for `gt_pca` objects
#'
#' For `gt_pca`, the following types of plots are available:
#' - `screeplot`: a plot of the eigenvalues of the principal components
#'    (currently it plots the singular value)
#' - `scores` a scatterplot of the scores of each individual on two principal
#' components (defined by `pc`)
#' - `loadings` a plot of loadings of all loci for a given component
#'    (chosen with `pc`)
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are not
#' customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' @param object an object of class `gt_pca`
#' @param type the type of plot (one of "screeplot", "scores" and "loadings")
#' @param k the principal components to be plotted: for scores, a pair of values
#'   e.g. c(1,2); for `loadings` either one or more values.
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @name autoplot_gt_pca
#' @export
#' @examples
#' library(ggplot2)
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' lobsters <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Remove monomorphic loci and impute
#' lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
#' lobsters <- gt_impute_simple(lobsters, method = "mode")
#'
#' # Create PCA object
#' pca <- gt_pca_partialSVD(lobsters)
#'
#' # Screeplot
#' autoplot(pca, type = "screeplot")
#'
#' # Scores plot
#' autoplot(pca, type = "scores")
#'
#' # Colour by population
#' autoplot(pca, type = "scores") + aes(colour = lobsters$population)
#'
#' # Scores plot of different components
#' autoplot(pca, type = "scores", k = c(1, 3)) +
#'   aes(colour = lobsters$population)
autoplot.gt_pca <- function(
    object,
    type = c("screeplot", "scores", "loadings"),
    k = NULL,
    ...) {
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (type == "screeplot") {
    plot(object, type = "screeplot")
  } else if (type == "scores") {
    if (is.null(k)) {
      k <- c(1, 2)
    }
    if (length(k) != 2) {
      stop("for 'scores' plots, 'pc' should be a pair of values, e.g. c(1,2)")
    }
    if ("square_frobenius" %in% names(object)) {
      tidy_pca <- tidy(object, "eigenvalues")
      k_percentage <- tidy_pca$percent[k]
      plot(object, type = "scores", scores = k) +
        ggplot2::labs(
          x = paste0(
            "PC",
            k[1],
            " (",
            round(k_percentage[1], 2),
            "%)"
          ),
          y = paste0("PC", k[2], " (", round(k_percentage[2], 2), "%)")
        )
    } else {
      plot(object, type = "scores", scores = k)
    }
  } else if (type == "loadings") {
    if (is.null(k)) {
      k <- 1
    }
    plot(object, type = "loadings", loadings = k)
  }
}
