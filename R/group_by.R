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
      c("gen_tbl", "grouped_gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")]) # nolint
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


#' @export
#' @noRd
dplyr_reconstruct.gen_tbl <- function(data, template) {
  out <- NextMethod()
  # if the genotypes are gone, drop the tbl_df class
  if (!"genotypes" %in% names(data)) {
    message("as genotypes were dropped, this is not longer a 'gen_tbl'")
    class(out) <- class(out)[-1]
  }
  out
}

#' @export
#' @noRd
dplyr_reconstruct.grouped_gen_tbl <- function(data, template) {
  out <- NextMethod()
  # if the genotypes are gone, drop the tbl_df class
  if (!"genotypes" %in% names(data)) {
    message("as genotypes were dropped, this is not longer a 'gen_tbl'")
    class(out) <- c("grouped_df", "tbl_df", "tbl", "data.frame")
  }
  out
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
