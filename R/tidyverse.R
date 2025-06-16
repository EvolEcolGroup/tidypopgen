#' Tidyverse methods for gt objects
#'
#' A filter method for `gen_tibble` objects
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#' @returns a `gen_tibble`
#' @export
#' @examples
#' test_gt <- example_gt("gen_tbl")
#' test_gt %>% filter(id %in% c("a", "c"))
#' test_gt <- example_gt("grouped_gen_tbl_sf")
#' test_gt %>% filter(id %in% c("a", "c"))
filter.gen_tbl <- function(..., deparse.level = 1) { # nolint
  # send it to the next (data.frame) method
  out <- NextMethod()
  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(out) <- obj_class
  }
  return(out)
}

#' An arrange method for `gen_tibble` objects
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#' @returns a `gen_tibble`
#' @export
#' @examples
#' test_gt <- example_gt("gen_tbl")
#' test_gt %>% arrange(id)
#' test_gt <- example_gt("grouped_gen_tbl_sf")
#' test_gt %>% arrange(id)
arrange.gen_tbl <- function(..., deparse.level = 1) { # nolint
  # send it to the next (data.frame) method
  out <- NextMethod()
  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(out) <- obj_class
  }
  return(out)
}


#' A mutate method for `gen_tibble` objects
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#' @returns a `gen_tibble`
#' @export
#' @examples
#' example_gt <- example_gt("gen_tbl")
#'
#' # Add a new column
#' example_gt %>% mutate(region = "East")
mutate.gen_tbl <- function(..., deparse.level = 1) { # nolint
  # send it to the next (data.frame) method
  out <- NextMethod()
  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(out) <- obj_class
  }
  return(out)
}

#' A dplyr_reconstruct method for `gen_tibble` objects
#'
#' @param data a gen_tibble and a data.frame or tibble
#' @param template a data.frame or tibble
#' @returns a `gen_tibble`
#' @export
dplyr_reconstruct.gen_tbl <- function(data, template) {
  # send it to the next (data.frame) method
  out <- NextMethod()
  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(out) <- obj_class
  }
  return(out)
}

#' A $ method for `gen_tibble` objects
#' @param data a gen_tibble
#' @param i column name
#' @param value a value to assign
#' @returns a `gen_tibble`
#' @export
#' @examples
#' example_gt <- example_gt("gen_tbl")
#'
#' # Add a new column
#' example_gt$region <- "East"
#'
#' example_gt
"$<-.gen_tbl" <- function(data, i, value) {
  # send it to the next (data.frame) method
  out <- NextMethod()
  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(out) <- obj_class
  }
  return(out)
}
