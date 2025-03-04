#' Autoplots for `gt_pcadapt` objects
#'
#' For `gt_pcadapt`, the following types of plots are available:
#' - `qq`: a qunatile-quantile plot of the p-values from `pcadapt`
#' (wrapping [bigsnpr::snp_qq()])
#' - `manhattan` a manhattan plot of the p-values from `pcadapt`
#' (wrapping [bigsnpr::snp_manhattan()])
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are
#' not customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' @param object an object of class `gt_pcadapt`
#' @param type the type of plot (one of "qq", and "mnahattan")
#' @param ... further arguments to be passed to [bigsnpr::snp_qq()] or
#' [bigsnpr::snp_manhattan()].
#' @returns a `ggplot2` object
#' @name autoplot_gt_pcadapt
#' @export
autoplot.gt_pcadapt <- function(object,
                                type = c("qq", "manhattan"),
                                ...) {
  type <- match.arg(type)
  if (type == "qq") {
    bigsnpr::snp_qq(object, ...)
  } else if (type == "manhattan") {
    bigsnpr::snp_manhattan(object,
      infos.chr = attr(object, "loci")$chromosome,
      infos.pos = attr(object, "loci")$position,
      ...
    )
  }
}
