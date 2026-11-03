#' Scale constructor using the distruct colours
#'
#' A wrapper around [ggplot2::scale_fill_manual()], using the distruct colours
#' from [`distruct_colours`].
#' @param guide guide function passed to [ggplot2::scale_fill_manual()].
#' Defaults to "none", set to "legend" if a legend is required.
#' @param ... further parameters to be passed to [ggplot2::scale_fill_manual()]
#' @returns a scale constructor to be used with ggplot
#' @export
#' @seealso [ggplot2::scale_fill_manual()] which this function wraps.
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
#' # Fill bars of individuals per population with the distruct colours
#' ggplot(lobsters, aes(x = population, fill = population)) +
#'   geom_bar() +
#'   scale_fill_distruct()
scale_fill_distruct <- function(guide = "none", ...) {
  ggplot2::scale_fill_manual(
    values = tidypopgen::distruct_colours,
    guide = guide,
    ...
  )
}
