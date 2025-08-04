# Note that we need to make sure that we deal with sf objects correctly, so
# that the sf methods don't override the gen_tbl

#' Tidyverse methods for gt objects
#'
#' A filter method for `gen_tibble` objects
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#' @returns a `gen_tibble`
#' @export
#' @examples
#' test_gt <- load_example_gt("gen_tbl")
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

#' A filter method for grouped `gen_tibble` objects
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#' @returns a grouped `gen_tibble`
#' @export
#' @examples
#' test_gt <- load_example_gt("grouped_gen_tbl")
#' test_gt %>% filter(id %in% c("a", "c"))
#' test_gt <- load_example_gt("grouped_gen_tbl_sf")
#' test_gt %>% filter(id %in% c("a", "c"))
filter.grouped_gen_tbl <- function(..., deparse.level = 1) { # nolint
  # send it to the next (data.frame) method
  out <- NextMethod()

  # prioritise "grouped_gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf", obj_class[!obj_class %in% c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf")]) # nolint
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
#' test_gt <- load_example_gt("gen_tbl")
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

#' An arrange method for grouped `gen_tibble` objects
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#' @returns a grouped `gen_tibble`
#' @export
#' @examples
#' test_gt <- load_example_gt("grouped_gen_tbl")
#' test_gt %>% arrange(id)
#' test_gt <- load_example_gt("grouped_gen_tbl_sf")
#' test_gt %>% arrange(id)
arrange.grouped_gen_tbl <- function(..., deparse.level = 1) { # nolint
  # send it to the next (data.frame) method
  out <- NextMethod()
  # prioritise "grouped_gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf", obj_class[!obj_class %in% c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf")]) # nolint
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
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Add a new column
#' example_gt %>% mutate(region = "East")
mutate.gen_tbl <- function(..., deparse.level = 1) { # nolint
  # send it to the next (data.frame) method
  out <- NextMethod()
  # the following check is a precaution, the message should not be triggered,
  # do not remove
  if (!"genotypes" %in% names(out)) {
    message("as genotypes were dropped, this is no longer a 'gen_tbl'")
    class(out) <- c("grouped_df", "tbl_df", "tbl", "data.frame")
    return(out)
  }
  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(out) <- obj_class
  }
  return(out)
}

#' A mutate method for grouped `gen_tibble` objects
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#' @returns a grouped `gen_tibble`
#' @export
#' @examples
#' test_gt <- load_example_gt("grouped_gen_tbl")
#' test_gt %>% mutate(region = "East")
#' test_gt <- load_example_gt("grouped_gen_tbl_sf")
#' test_gt %>% mutate(region = "East")
mutate.grouped_gen_tbl <- function(..., deparse.level = 1) { # nolint
  # send it to the next (data.frame) method
  out <- NextMethod()
  # the following check is a precaution, the message should not be triggered,
  # do not remove
  if (!"genotypes" %in% names(out)) {
    message("as genotypes were dropped, this is no longer a 'gen_tbl'")
    class(out) <- c("grouped_df", "tbl_df", "tbl", "data.frame")
    return(out)
  }
  # prioritise "grouped_gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf", obj_class[!obj_class %in% c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf")]) # nolint
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
#' @noRd
dplyr_reconstruct.gen_tbl <- function(data, template) {
  # send it to the next (data.frame) method
  out <- NextMethod()
  # the following check is a precaution, the message should not be triggered,
  # do not remove
  if (!"genotypes" %in% names(data)) {
    message("as genotypes were dropped, this is no longer a 'gen_tbl'")
    class(out) <- c("tbl_df", "tbl", "data.frame")
    return(out)
  }
  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(out) <- obj_class
  }
  return(out)
}

#' A dplyr_reconstruct method for grouped `gen_tibble` objects
#'
#' @param data a gen_tibble and a data.frame or tibble
#' @param template a data.frame or tibble
#' @returns a grouped `gen_tibble`
#' @export
#' @noRd
dplyr_reconstruct.grouped_gen_tbl <- function(data, template) {
  out <- NextMethod()
  # the following check is a precaution, the message should not be triggered,
  # do not remove
  if (!"genotypes" %in% names(data)) {
    message("as genotypes were dropped, this is no longer a 'gen_tbl'")
    class(out) <- c("grouped_df", "tbl_df", "tbl", "data.frame")
    return(out)
  }
  # prioritise "grouped_gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf", obj_class[!obj_class %in% c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf")]) # nolint
    class(out) <- obj_class
  }
  return(out)
}

#' A $ method for `gen_tibble` objects
#' @param x a gen_tibble
#' @param i column name
#' @param value a value to assign
#' @returns a `gen_tibble`
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Add a new column
#' example_gt$region <- "East"
#'
#' example_gt
"$<-.gen_tbl" <- function(x, i, value) {
  if (is.null(value) && inherits(x[[i]], "sfc") && ((is.character(i) && # nolint start
    i == attr(x, "sf_column")) || (is.integer(i) && names(x)[i] ==
    attr(x, "sf_column")))) { # nolint end
    sf::st_set_geometry(x, NULL)
  } else {
    x[[i]] <- value
    x
  }
  return(x)
}

#' @export
#' @noRd
group_by.gen_tbl <- function(
    .data,
    ...,
    .add = FALSE,
    .drop = group_by_drop_default(.data)) {
  out <- NextMethod()

  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(out)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf", obj_class[!obj_class %in% c("grouped_df", "gen_tbl", "sf")]) # nolint
    class(out) <- obj_class
  } else {
    class(out) <- c("grouped_gen_tbl", "grouped_df", "gen_tbl", class(out)[-1])
  }
  return(out)
}

#' @export
#' @noRd
ungroup.grouped_gen_tbl <- function(x, ...) {
  out <- NextMethod(...)
  class(out) <- c("gen_tbl", class(out))
  return(out) # nolint
}

# drop the `gen_tbl` class if the `genotype` column is subsetted out
#' @export
#' @noRd
"[.gen_tbl" <- function(x, i, j, ...) {
  x <- NextMethod()
  if (!"genotypes" %in% names(x)) {
    class(x) <- class(x)[!class(x) %in% "gen_tbl"]
  }
  x
}

# drop the `gen_tbl` class if the `genotype` column is subsetted out
#' @export
#' @noRd
"[.grouped_gen_tbl" <- function(x, i, j, ...) {
  original_class <- class(x)
  x <- NextMethod()
  class(x) <- original_class
  if (!"genotypes" %in% names(x)) {
    class(x) <- class(x)[!class(x) %in% c("grouped_gen_tbl", "gen_tbl")]
  }
  x
}

# drop the `gen_tbl` class if the `genotype` column is subsetted out
#' @export
#' @noRd
"[[<-.gen_tbl" <- function(x, i, value) {
  x <- NextMethod()
  if (!"genotypes" %in% names(x)) {
    class(x) <- class(x)[!class(x) %in% "gen_tbl"]
  }
  obj_class <- class(x)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(x) <- obj_class
  }
  x
}

# drop the `gen_tbl` class if the `genotype` column is subsetted out
#' @export
#' @noRd
"[[<-.grouped_gen_tbl" <- function(x, i, value) {
  original_class <- class(x)
  x <- NextMethod()
  class(x) <- original_class
  if (!"genotypes" %in% names(x)) {
    class(x) <- class(x)[!class(x) %in% c("grouped_gen_tbl", "gen_tbl")]
  }
  # prioritise "grouped_gen_tbl" class over "sf"
  obj_class <- class(x)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf", obj_class[!obj_class %in% c("grouped_gen_tbl", "grouped_df", "gen_tbl", "sf")]) # nolint
    class(x) <- obj_class
  }
  x
}


#' @export
#' @noRd
dplyr_row_slice.grouped_gen_tbl <- function(data, i, ...) {
  original_class <- class(data)
  x <- NextMethod()
  class(x) <- original_class
  x
}

#' @export
#' @noRd
dplyr_col_modify.grouped_gen_tbl <- function(data, cols) {
  original_class <- class(data)
  x <- NextMethod()
  class(x) <- original_class
  x
}

# nolint start
# #' @export
# dplyr_row_slice.gen_tbl<-function(data, i, ...){
#   NextMethod()
#
# }

# #' @export
# dplyr_col_modify.gen_tbl<-function(data, cols){
#   NextMethod()
#
# }
# nolint end
