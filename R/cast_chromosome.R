cast_chromosome_to_factor <- function(x){
  if(is.numeric(x)){
    x <- as.factor(x)
  }
  if(is.character(x)){
    # prepare for transforming to a factor
    digits <- x %>% stringr::str_match_all("\\d+") %>% sapply(toString)
    non_digits <- x[!stringr::str_detect(x, "^[0-9]+$")]

    # transform digits to numeric first (so as.factor orders them correctly)
    digits <- as.numeric(digits)

    # then transform to a factor
    digits_fac <- as.factor(digits)

    # find NA positions (those which were non-digits originally)
    na_pos <- which(is.na(digits))

    # fill these with their original values
    digits[na_pos] <- non_digits

    # then use levels to construct a factor for all values
    x <- factor(digits,
                          levels = c(levels(digits_fac), unique(non_digits)))

  }
  return(x)
}

cast_chromosome_to_int <- function(chromosome) {
  if(length(stringr::str_detect(chromosome, "^[0-9]+$") >1)){
    digits <- chromosome %>% stringr::str_match_all("\\d+") %>% sapply(toString)
    non_digits <- chromosome[!stringr::str_detect(chromosome, "^[0-9]+$")]
    non_digits <- droplevels(non_digits)

    # transform digits to integer first (so orders them with correct value,
    # not with the levels of the factor)
    digits <- as.integer(digits)

    # find highest value of digits
    max_digit <- max(digits, na.rm=TRUE)
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
