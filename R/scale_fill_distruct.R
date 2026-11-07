#' Scale constructor using the distruct colours
#'
#' A discrete fill scale using the distruct colours from
#' [`distruct_colours`], in order. With more groups than colours (56), the
#' colours are recycled, with a warning, so some groups share a colour.
#' @param guide guide function passed to [ggplot2::discrete_scale()].
#' Defaults to "none", set to "legend" if a legend is required.
#' @param ... further parameters to be passed to [ggplot2::discrete_scale()]
#'   (e.g. `name`, `labels` or `breaks`)
#' @returns a scale constructor to be used with ggplot
#' @export
#' @seealso [ggplot2::scale_fill_manual()], which this function used to wrap
#'   (it can not recycle colours).
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
  ggplot2::discrete_scale(
    aesthetics = "fill",
    palette = distruct_palette,
    guide = guide,
    # the default of ggplot2::scale_fill_manual()
    na.value = "grey50",
    ...
  )
}

# The first n distruct colours, recycled (with a warning) if there are more
# than the number of colours
distruct_palette <- function(n) {
  colours <- tidypopgen::distruct_colours
  if (n > length(colours)) {
    warning(
      "there are ", n, " groups, but only ", length(colours),
      " distruct colours: colours are recycled, so some groups share a colour"
    )
  }
  rep_len(colours, n)
}
