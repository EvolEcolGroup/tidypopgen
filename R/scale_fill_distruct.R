#' Scale constructor using the distruct colours
#'
#' A wrapper around [ggplot2::scale_fill_manual()], using the distruct colours
#' from [`distruct_colours`].
#' @param guide guide function passed to [ggplot2::scale_fill_manual()].
#' Defaults to "none", set to "legend" if a legend is required.
#' @param ... further parameters to be passed to [ggplot2::scale_fill_manual()]
#' @returns a scale constructor to be used with ggplot
#' @export

scale_fill_distruct <- function(guide = "none", ...) {
  ggplot2::scale_fill_manual(
    values = tidypopgen::distruct_colours,
    guide = guide, ...
  )
}
