# replace
# forcats::fct_inorder
fct_inorder_base <- function(f, ordered = FALSE) {
  if (!is.factor(f)) stop("Input must be a factor")
  if (!is.logical(ordered) || length(ordered) != 1) stop("ordered must be a single logical value")

  idx <- match(levels(f), f, nomatch = 0) # Find first occurrences of levels
  idx <- unique(idx[idx > 0]) # Remove missing matches
  idx <- union(idx, seq_along(levels(f))) # Ensure all levels are included

  factor(f, levels = levels(f)[idx], ordered = ordered)
}

# Example usage
x <- factor(c("b", "a", "c", "a", "b", "c", "b", "a"))
testthat::expect_identical(forcats::fct_inorder(x), fct_inorder_base(x))


########
# replace stringr::str_replace
str_replace_base <- function(string, pattern, replacement) {
  sub(pattern, replacement, string)
}

# Example usage:
x <- "I love cats and cats are cute"
testthat::expect_identical(stringr::str_replace(x, "cats", "dogs"),
                           str_replace_base(x, "cats", "dogs"))

