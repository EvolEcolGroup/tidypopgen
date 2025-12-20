# utility internal functions


# simple function to extract the extension of a file
# see tools::file_ext
file_ext <- function(x) {
  utils::tail(unlist(strsplit(x, ".", fixed = TRUE)), n = 1)
}


#######
# Convenience functions that are not exported by `bigstatsr`
CutBySize <- function(m, block.size, nb = ceiling(m / block.size)) { # nolint
  bigparallelr::split_len(m, nb_split = nb)
}

seq2 <- function(lims) {
  seq(lims[1], lims[2])
}
#######


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

  # idx <- match(levels(f), f, nomatch = 0) # Find first occurrences of levels #nolint start
  # idx <- unique(idx[idx > 0]) # Remove missing matches
  # idx <- union(idx, seq_along(levels(f))) # Ensure all levels are included #nolint end

  idx <- as.integer(f)[!duplicated(f)]
  idx <- union(idx[!is.na(idx)], seq_along(levels(f)))

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
