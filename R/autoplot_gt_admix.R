#' Autoplots for `gt_admix` objects
#'
#' For `gt_admix`, the following types of plots are available:
#' - `cv`: the cross-validation error for each value of `k`
#' - `barplot` a standard barplot of the admixture proportions
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are
#' not customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' This autoplot will automatically rearrange individuals according to
#' their id and any grouping variables if an associated 'data' gen_tibble is provided.
#' To avoid any automatic re-sorting of individuals, set `arrange_by_group` and
#' `arrange_by_indiv` to FALSE. See `autoplot.q_matrix` for further details.
#'
#' @param object an object of class `gt_admixture`
#' @param type the type of plot (one of "cv", and "boxplot")
#' @param k the value of `k` to plot (for `barplot` type only)
#' param repeat the repeat to plot (for `barplot` type only)
#' @param run the run to plot (for `barplot` type only)
#' @param ... additional arguments to be passed to autoplot method
#' for q_matrices [autoplot_q_matrix()], used when type is `barplot`.
#' @returns a `ggplot2` object
#' @name autoplot_gt_admix
#' @export
autoplot.gt_admix <- function(object,
                            type=c("cv", "barplot"),
                            k = NULL,
                            run = NULL,
                            ...)
{
  type <- match.arg(type)
  if (type== "cv") {
    if (is.null(object$cv)){
      stop("No cross validation error available")
    }

    if(object$algorithm == "SNMF"){
      ggplot2::ggplot(data.frame(k=object$k, cv=object$cv), ggplot2::aes(x=.data$k, y=.data$cv)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(x="k", y="Cross-Entropy")
    } else {
      ggplot2::ggplot(data.frame(k=object$k, cv=object$cv), ggplot2::aes(x=.data$k, y=.data$cv)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(x="k", y="Cross validation error")
    }

  } else if (type == "barplot") {
    # check that k is specified
    if (is.null(k)){
      stop("You must specify a value for k")
    }
    # check that run is specified
    if (is.null(run)){
      stop("You must specify a value for repeat")
    }
    # get the Q matrix for the specified k and repeat
    Q <- get_q_matrix(object, k = k, run = run)
    autoplot(Q, ...)
  }
}
