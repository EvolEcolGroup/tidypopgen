#' Add an simple feature geometry to a `gen_tibble`
#'
#' `gt_add_sf` adds an active sf geometry column to a `gen_tibble` object.  The
#' resulting `gen_tbl` inherits from `sf` and can be used with functions from
#' the `sf` package. It is possible to either create a [`sf::sfc`] geometry
#' column from coordinates, or to provide an existing geometry column (which
#' will then become the active geometry for `sf`).
#'
#' @param x a [`gen_tibble`] object
#' @param coords a vector of length 2, giving the names of the x and y columns
#'   in `x` (i.e. the coordinates, e.g. longitude and latitude). If `coords` is
#'   not provided, the geometry column must be provided.
#' @param crs the coordinate reference system of the coordinates. If this is not
#'   set, it will be set to the default value of `sf::st_crs(4326)`.
#' @param sfc_column the name of an [`sf::sfc`] column to be used as the
#'   geometry
#' @return a [`gen_tibble`] object with an additional geometry column (and thus
#'   belonging also to `sf` class).
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Add some coordinates
#' example_gt <- example_gt %>% mutate(
#'   longitude = c(0, 0, 2, 2, 0, 2, 2),
#'   latitude = c(51, 51, 49, 49, 51, 41, 41)
#' )
#'
#' # Convert lat and long to sf:
#' example_gt <- gt_add_sf(x = example_gt, coords = c("longitude", "latitude"))
#'
#' # Check class
#' class(example_gt)
gt_add_sf <- function(x, coords = NULL, crs = NULL, sfc_column = NULL) {
  # check that x is a gen_tibble
  if (!inherits(x, "gen_tbl")) {
    stop("x must be a gen_tibble")
  }
  # check that x does not already inherit from sf
  if (inherits(x, "sf")) {
    stop("x already inherits from sf")
  }
  # check that coords and geometry are not both provided
  if (!is.null(coords) && !is.null(sfc_column)) {
    stop("You must provide either coords or sfc_column, not both")
  }
  # check that at least one of coords or geometry is provided
  if (is.null(coords) && is.null(sfc_column)) {
    stop("You must provide either coords or sfc_column")
  }
  # check that crs is a valid crs
  if (!is.null(crs)) {
    if (!inherits(crs, "crs")) {
      stop("crs must be a valid crs object")
    }
  } else {
    # set crs to default value
    crs <- sf::st_crs(4326)
  }
  # check that coords includes two columns from x
  if (!is.null(coords)) {
    if (length(coords) != 2) {
      stop("coords must be a vector of length 2")
    }
    if (!all(coords %in% colnames(x))) {
      stop("coords must be a vector of length 2, with names of columns in x")
    }
    # create geometry column from coords
    geometry_sf <- sf::st_as_sf(
      x %>% select(dplyr::all_of(coords)),
      coords = coords,
      crs = crs,
      remove = TRUE
    )
    # add geometry to x
    if ("geometry" %in% colnames(x)) {
      warning(
        "x already has a column named 'geometry'; this will ",
        "be overwritten"
      )
    }
    x$geometry <- sf::st_geometry(geometry_sf)
    sfc_column <- "geometry"
  }

  # check that sfc_column is a valid geometry
  if (!is.null(sfc_column)) {
    # check that the geometry column exists in x
    if (!sfc_column %in% colnames(x)) {
      stop(paste0("sfc_column '", sfc_column, "' does not exist in x"))
    }
    # and that it is an sfc object
    if (!inherits(x[[sfc_column]], "sfc")) {
      stop("sfc_column must be a valid sfc object")
    }
  }

  # add attribute sf_column
  attr(x, "sf_column") <- sfc_column
  # add sf class
  new_class_def <- c("gen_tbl", "sf", "tbl_df", "tbl", "data.frame")
  if (inherits(x, "grouped_gen_tbl")) {
    # if x is a grouped gen_tbl, we need to add the group classes in front
    new_class_def <- c("grouped_gen_tbl", "grouped_df", new_class_def)
  }
  class(x) <- new_class_def
  return(x)
}
