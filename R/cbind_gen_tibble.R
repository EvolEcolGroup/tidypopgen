#' Combine a gen_tibble to a data.frame or tibble by column
#'
#' A [cbind()] method to merge `gen_tibble` objects with data.frames and normal
#' tibbles. Whilst this works, it is not ideal as it does not check the order of
#' the tables, and we suggest that you use [`dplyr::left_join()`] instead. Note
#' that `cbind` will not combine two `gen_tibbles` (i.e. it will NOT combine
#' markers for the same individuals)
#'
#' As in [dplyr::bind_cols()] (which `dplyr` uses to `cbind()` grouped
#' tibbles), duplicated column names are made unique (e.g. `population...2`
#' and `population...4`), with a message. If a grouping variable is renamed,
#' the result is no longer grouped. The columns `id` and
#' `genotypes` of the `gen_tibble` can not be duplicated, as renaming them
#' would break it.
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#'   See \code{\link[base]{cbind}} for details.
#' @returns a `gen_tibble`
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Create a dataframe to combine with the gen_tibble
#' df <- data.frame(region = c("A", "A", "B", "B", "A", "B", "B"))
#'
#' # Combine the gen_tibble with the dataframe
#' example_gt <- cbind(example_gt, df)
cbind.gen_tbl <- function(..., deparse.level = 1) { # nolint
  args <- list(...)
  which_gen_tbl <- which(sapply(args, inherits, "gen_tbl"))[1]
  class_orig <- class(args[[which_gen_tbl]])
  stop_if_cbind_reserved(args, which_gen_tbl)

  # send it to the data.frame method
  out <- cbind.data.frame(..., deparse.level = deparse.level)
  # make duplicated names unique, with a message, as dplyr::bind_cols() does
  # (a grouped gen_tibble has already been repaired by dplyr)
  names(out) <- vctrs::vec_as_names(names(out), repair = "unique")
  # sort out the lost class information
  class(out) <- class_orig
  # prioritise "gen_tbl" class over "sf"
  if ("sf" %in% class_orig) {
    class_orig <-
      c("gen_tbl", "sf", class_orig[!class_orig %in% c("gen_tbl", "sf")])
    class(out) <- class_orig
  }
  return(out)
}

#' @rdname cbind.gen_tbl
#' @export
cbind.grouped_gen_tbl <- function(..., deparse.level = 1) { # nolint
  args <- list(...)
  stop_if_cbind_reserved(args, which(sapply(args, inherits, "gen_tbl"))[1])
  # as dplyr does for grouped data.frames (cbind.grouped_df())
  dplyr::bind_cols(...)
}

# stop if the columns of the arguments of cbind() other than the gen_tibble
# (the names of data.frames and matrices, or the argument names of vectors)
# include the id or the genotypes, as renaming those would break the
# gen_tibble
stop_if_cbind_reserved <- function(args, which_gen_tbl) {
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }
  new_cols <- unlist(lapply(seq_along(args)[-which_gen_tbl], function(i) {
    if (is.data.frame(args[[i]])) {
      names(args[[i]])
    } else if (!is.null(colnames(args[[i]]))) {
      colnames(args[[i]])
    } else {
      arg_names[i]
    }
  }))
  reserved <- intersect(c("id", "genotypes"), new_cols)
  if (length(reserved) > 0) {
    stop(
      "the gen_tibble already has a column ",
      paste0("'", reserved, "'", collapse = " and "),
      ", which can not be duplicated: rename it in the other object, ",
      "or use dplyr::left_join()",
      call. = FALSE
    )
  }
  invisible(NULL)
}
