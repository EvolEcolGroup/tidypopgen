# utility internal functions

# simple function to extract the extension of a file
# this avoids having to import tools just for tools::file_ext
file_ext <- function(x) {
  utils::tail(unlist(strsplit(x, ".", fixed = TRUE)), n = 1)
}

# a simple tidier for dist matrices
# tidy.dist is deprecated, and often we have a full matrix rather
# than a dist object
tidy_dist_matrix <- function(mat) {
  if (!inherits(mat, "matrix")) {
    stop("mat should be a matrix")
  }
  xy <- t(utils::combn(colnames(mat), 2))
  colnames(xy) <- c("item1", "item2")
  xy %>%
    tibble::as_tibble() %>%
    dplyr::mutate(value = mat[xy])
}

# stop if not diploid
stopifnot_diploid <- function(x) {
  if (inherits(x, "gen_tbl")) {
    x <- x$genotypes
  }
  if (attr(x, "ploidy") == -2) {
    if (min(attr(x, "bigsnp")$fam$ploidy[vctrs::vec_data(x)]) != 2) {
      stop("this function only works on diploid data")
    }
  }
}

stopifnot_dip_pseudo <- function(x) {
  if (inherits(x, "gen_tbl")) {
    x <- x$genotypes
  }
  if (abs(attr(x, "ploidy")) != 2) {
    stop("this function only works on diploid data")
  }
}


is_diploid_only <- function(x) {
  if (inherits(x, "gen_tbl")) {
    (attr(x$genotypes, "ploidy") == 2)
  } else {
    (attr(x, "ploidy") == 2)
  }
}

is_pseudohaploid <- function(x) {
  if (inherits(x, "gen_tbl")) {
    (attr(x$genotypes, "ploidy") == -2)
  } else {
    (attr(x, "ploidy") == -2)
  }
}


#' Reorder factor levels by first appearance
#'
#' This function reorders the levels of a factor so that they appear in the
#' order they first appear in the data.  This is a drop-in replacement for
#' `forcats::fct_inorder()`
#'
#' @param f A factor or character vector
#' @param ordered Logical, should the resulting factor be ordered?
#' @return A factor with levels in the order they first appear
#' @keywords internal
#' @noRd
# nolint start
# x <- factor(c("b", "a", "c", "a", "b", "c", "b", "a"))
# fct_inorder_base(x)
# fct_inorder_base(x, ordered = TRUE)
# nolint end

fct_inorder_base <- function(f, ordered = FALSE) {
  if (is.character(f)) {
    f <- factor(f)
  }
  if (!is.logical(ordered) || length(ordered) != 1) {
    stop("ordered must be a single logical value")
  }

  idx <- match(levels(f), f, nomatch = 0) # Find first occurrences of levels
  idx <- unique(idx[idx > 0]) # Remove missing matches
  idx <- union(idx, seq_along(levels(f))) # Ensure all levels are included

  factor(f, levels = levels(f)[idx], ordered = ordered)
}

#' Replace matches with new text
#'
#' This function is a drop-in replacement for `stringr::str_replace()`, which
#' replaces the first match
#'
#' @param string A character input vector
#' @param pattern A regular expression
#' @param replacement A character vector
#' @return A character vector with the first match replaced
#' @keywords internal
#' @noRd
# nolint start
# x <- "I love cats and cats are cute"
# str_replace_base(x, "cats", "dogs")
# nolint end

str_replace_base <- function(string, pattern, replacement) {
  sub(pattern, replacement, string)
}
