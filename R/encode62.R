#' Function to encode a pair of integers into a base62 string.
#'
#' This function is used to encode the chromosome and the coordinates of a SNP
#' into a compact string representation using base62 encoding. The first integer
#' represents the chromosome number, and the second integer represents the
#' position on that chromosome. The encoded string can be used as a unique
#' identifier for the SNP.
#'
#' The two numbers are combined by exploiting the fact that the the number of
#' chromosome is usually small and known in advance (e.g., 24 for humans). This
#' allows the chromosome number to be encoded using a fixed number of base-62
#' digits, while the position can use a variable number of digits. This ensures
#' that the combined encoding is unambiguous and compact. Specifically, both
#' numbers are converted to base-62 using only ASCII characters. The position is
#' encoded first, using as many base-62 digits as needed. The first number is
#' then encoded using a fixed number of base-62 digits (chosen in advance to be
#' large enough to represent its maximum possible value) and left-padded with
#' zeros if necessary. The final string is formed by concatenating these two
#' encodings, with the fixed-width encoding of the first number placed at the
#' end. During decoding, the last fixed number of characters are read back as
#' the chromosome number, and the remaining prefix is decoded as the position,
#' allowing the original pair to be recovered unambiguously without any
#' separator characters.
#'
#' @param chr Chromosomes, either character/factor, which will be cast to
#'   integer, or integer representing the chromosome number.
#' @param pos Integer representing the position on the chromosome.
#' @param max_chr Integer representing the maximum chromosome number expected.
#'   This is used to determine the fixed width for encoding the chromosome
#'   number. If left NULL, the maximum value of `chr` will be used.
#' @returns A vector of base62 encoded strings representing the combined
#'   chromosomes and positions, with the maximum chromosome number stored as an
#'   attribute named `encode62_max_chr`.
#' @export
#' @examples
#' encoded_coords <- encode62(c(1, 10, 260), c(1, 1000, 1000000))
#' print(encoded_coords)
#' decoded_coords <- decode62(encoded_coords)
#' print(decoded_coords)
encode62 <- function(chr, pos, max_chr = NULL) {
  # check that chromosome and position are integer vectors of the same length
  if (length(chr) != length(pos)) {
    stop("chr and pos must be integer vectors of the same length")
  }
  # if chr is a string, turn it into a factor
  if (is.character(chr)) {
    chr <- as.factor(chr)
  }
  # convert chr to integer
  chr <- as.integer(chr)
  # convert pos to integer
  pos <- as.integer(pos)

  if (is.null(max_chr)) {
    max_chr <- max(chr)
  }
  # check number of digits needed for chr
  num_digits_chr <- ceiling(log(max_chr + 1, base = 62))
  # use the cpp encoding function
  encode_vec <- encode_pair_vec_cpp(chr, pos, num_digits_chr)
  attr(encode_vec, "recode62_max_chr") <- max_chr
  return(encode_vec)
}

#' Function to decode a base62 encoded string back into a pair of integers.
#'
#' This function reverses the encoding performed by `encode_62`, taking a base62
#' encoded string and extracting the original chromosome number and position. It
#' uses the fixed width determined by the maximum chromosome number to correctly
#' separate the two components of the encoding.
#' @param encoded_str A base62 encoded string representing the combined
#'   chromosome and position
#' @param max_chr Integer representing the maximum chromosome number expected.
#'   This is used to determine the fixed width for decoding the chromosome
#'   number. If left NULL, an error will be raised, unless encoded_str has an
#'   attribute `recode62_max_chr`, in which case that value will be used.
#' @returns A maxtrix with two columns: `chr` and `pos`, representing the
#'   decoded chromosome numbers and positions.
#' @export
#' @examples
#' encoded_coords <- encode62(c(1, 10, 260), c(1, 1000, 1000000))
#' print(encoded_coords)
#' decoded_coords <- decode62(encoded_coords)
#' print(decoded_coords)
decode62 <- function(encoded_str, max_chr = NULL) {
  # if max_chr is NULL, try to get it from the attribute
  if (is.null(max_chr)) {
    # check that the attribute exists, otherwise throw an error
    if (is.null(attr(encoded_str, "recode62_max_chr"))) {
      stop(
        "max_chr is NULL and encoded_str does not have ",
        "the attribute 'recode62_max_chr'"
      )
    }
    max_chr <- attr(encoded_str, "recode62_max_chr")
  }
  # estimate the number of digits needed for chr
  num_digits_chr <- ceiling(log(max_chr + 1, base = 62))
  # use the cpp decoding function
  decode_pair_vec_cpp(encoded_str, num_digits_chr)
}
