#' A function to transform all chromosome name inputs to a factor
#'
#' This function will order chromosome names correctly when they are numeric or
#' character vectors. It will handle mixed character vectors (e.g. "1", "2",
#' "X", "Y") by ordering numeric values first, then non-numeric values.
#'
#' @param x a vector of chromosome names (numeric or character)
#' @returns a factor with chromosome names ordered correctly
#' @keywords internal
#' @noRd
cast_chromosome_to_factor <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    x <- as.factor(x)
  } else if (is.character(x)) {
    # prepare for transforming to a factor
    digits <- x %>%
      base_str_match_all("\\d+") %>%
      sapply(toString)
    non_digits <- x[!grepl("^[0-9]+$", x)]

    # transform digits to numeric first (so as.factor orders them correctly)
    digits <- as.numeric(digits)

    # if all values are prefixed/suffixed with characters
    if (length(digits) == length(non_digits)) {
      # transform to a factor
      x <- as.factor(digits)

      # otherwise, if mixed character vector (e.g. 1, 2, "X", "Y")
    } else {
      # then transform to a factor
      digits_fac <- as.factor(digits)

      # find NA positions (those which were non-digits originally)
      na_pos <- which(is.na(digits))

      # fill these with their original values
      digits[na_pos] <- non_digits

      # then use levels to construct a factor for all values
      x <- factor(digits,
        levels = c(levels(digits_fac), unique(non_digits))
      )
    }
  }
  return(x)
}


#' A function to transform all chromosome name inputs to integers
#'
#' @param chromosome the chromosome column of a gen_tibble
#' @returns an integer vector with chromosome names converted to integers
#' @keywords internal
#' @noRd
cast_chromosome_to_int <- function(chromosome) {
  if (any(!grepl("^[0-9]+$", chromosome))) {
    digits <- chromosome %>%
      base_str_match_all("\\d+") %>%
      sapply(toString)
    non_digits <- chromosome[!grepl("^[0-9]+$", chromosome)]
    non_digits <- droplevels(non_digits)

    # transform digits to integer first (so orders them with correct value,
    # not with the levels of the factor)
    digits <- as.integer(digits)

    # find highest value of digits
    max_digit <- max(digits, na.rm = TRUE)
    # give non_digits a number higher than max_digit
    non_digit_values <- seq(max_digit + 1, max_digit + length(non_digits))

    # assign levels of non_digits to these new values
    levels(non_digits) <- non_digit_values

    # find NA positions (those which were non-digits originally)
    na_pos <- which(is.na(digits))

    # fill these with their original values
    digits[na_pos] <- as.integer(as.character(non_digits))

    chromosome <- digits
  } else {
    chromosome <- as.character(chromosome)
    chromosome <- as.integer(chromosome)
  }
  return(chromosome)
}
