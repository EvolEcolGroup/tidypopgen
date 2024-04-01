#' Autoplots for `gt_dapc` objects
#'
#' For `gt_dapc`, the following types of plots are available:
#' - `screeplot`: a plot of the eigenvalues of the discriminant axes
#' - `scores` a scatterplot of the scores of each individual on two discriminant
#' axes (defined by `ld`)
#' - `loadings` a plot of loadings of all loci for a discriminant axis (chosen with `ld`)
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are
#' not customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' @param object an object of class `gt_dapc`
#' @param type the type of plot (one of "screeplot", "scores" and "loadings")
#' @param ld the principal components to be plotted: for scores, a pair of values
#' e.g. c(1,2); for `loadings` either one or more values.
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @rdname autoplot_gt_pca
#' @export
autoplot.gt_dapc <- function(object,
                             type=c("screeplot", "scores","loadings"),
                             ld = NULL, ...)
{
  rlang::check_dots_empty()
  type <- match.arg(type)
  #  stop("autoplot for gt_dapc not avaialble yet")
  if (type== "screeplot") {
    tidy(object, matrix="eigenvalues") %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$LD,y=.data$eigenvalue)) +
      ggplot2::geom_point()+
      ggplot2::geom_line()
  } else if (type == "scores"){
    if (is.null(ld)){
      ld <- c(1,2)
    }
    if (length(ld)!=2){
      stop("for 'scores' plots, 'ld' should be a pair of values, e.g. c(1,2)")
    }
    tibble(cluster=object$grp) %>%
      mutate(LDa = object$ind.coord[,ld[1]],
             LDb = object$ind.coord[,ld[2]]) %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$LDa, y=.data$LDb,
                                   colour = .data$cluster))+
      ggplot2::geom_point()+
      ggplot2::stat_ellipse()+
      ggplot2::labs(x=paste0("LD",ld[1]), y=paste0("LD",ld[2]))
  } else if (type == "loadings"){
    if (is.null(ld)){
      ld <- 1
    }
    plot(object, type = "loadings", loadings = ld)
  }
}
