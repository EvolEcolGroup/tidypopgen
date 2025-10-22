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
  if (any(is.na(x))) {
    stop("NA values are not allowed in chromosome names")
  }
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

    # if all values are prefixed/suffixed with some characters
    if (length(which(!is.na(digits))) == length(non_digits)) {
      # transform to a factor
      dig_fac <- as.factor(digits)

      res <- extract_parts(x)

      # find unique entries for res - should be same length as dig_fac
      unique_res <- unique(res)

      # reorder res to match levels of dig_fac
      unique_res <- unique_res[match(levels(dig_fac), unique_res$digits), ]

      # reconstruct levels
      concat <- unique_res %>% tidyr::unite(concat,
        c(prefix, digits, suffix), # nolint
        sep = "",
        na.rm = TRUE
      )
      x <- factor(x, levels = concat$concat)

      # otherwise, if mixed character vector
      # (e.g. 1, 2, "X", "Y" or "chr1", "chr2", "chrX", "chrY")
    } else {
      # then transform to a factor
      digits_fac <- as.factor(digits)

      # list of non-digit unique entries
      no_digits <- x[!grepl("[0-9]", x)]

      res <- extract_parts(x)
      # find unique entries for res - should be same length as dig_fac
      unique_res <- unique(res)

      # reorder res to match levels of dig_fac
      unique_res <- unique_res[match(levels(digits_fac), unique_res$digits), ]

      # reconstruct levels
      concat <- unique_res %>% tidyr::unite(concat,
        c(prefix, digits, suffix), # nolint
        sep = "",
        na.rm = TRUE
      )

      # find NA positions (those which were non-digits originally)
      na_pos <- which(is.na(digits))

      # fill these with their original values
      digits[na_pos] <- no_digits

      no_digits <- as.factor(no_digits)

      # then use levels to construct a factor for all values
      x <- factor(x, levels = c(concat$concat, levels(no_digits)))
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

    # transform digits to integer first (so orders them with correct value,
    # not with the levels of the factor)
    digits <- as.integer(digits)

    # list of non-digit unique entries
    no_digits <- chromosome[!grepl("[0-9]", chromosome)]
    no_digits <- droplevels(no_digits)

    # if there are no numeric values in chromosome
    if (all(is.na(digits))) {
      max_digit <- 0
    } else {
      # find highest value of digits
      max_digit <- max(digits, na.rm = TRUE)
    }

    # give non_digits a number higher than max_digit
    non_digit_values <- seq(max_digit + 1, max_digit + length(no_digits))

    # assign levels of non_digits to these new values
    levels(no_digits) <- non_digit_values

    # find NA positions (those which were non-digits originally)
    na_pos <- which(is.na(digits))

    # fill these with their original values
    digits[na_pos] <- as.integer(as.character(no_digits))

    chromosome <- digits
  } else {
    chromosome <- as.character(chromosome)
    chromosome <- as.integer(chromosome)
  }
  return(chromosome)
}


extract_parts <- function(x) {
  # Regex: optional prefix (letters/underscores),
  # digits, optional suffix (letters/underscores)
  pattern <- "^([A-Za-z_]+)?(\\d+)([A-Za-z_]+)?$"

  # Match & extract pattern
  m <- regexec(pattern, x)
  parts <- regmatches(x, m)

  # Build tidy rows
  df_list <- lapply(parts, function(m) {
    if (length(m) == 4) {
      prefix <- ifelse(m[2] != "", m[2], NA)
      digits <- ifelse(m[3] != "", m[3], NA)
      suffix <- ifelse(m[4] != "", m[4], NA)

      data.frame(
        digits = as.numeric(digits),
        prefix = prefix,
        suffix = suffix,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        digits = NA_real_,
        prefix = NA_character_,
        suffix = NA_character_
      )
    }
  })

  df <- do.call(rbind, df_list)
  # keep rows with at least one non-NA part
  df <- df[!(is.na(df$digits) & is.na(df$prefix) & is.na(df$suffix)), , drop = FALSE] #nolint
  rownames(df) <- NULL
  df
}
