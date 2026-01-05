#' Function to encode a pair of integers into a base62 string.
#'
#' This function is used to encode the chromosome and the coordinates of a SNP
#' into a compact string representation using base62 encoding. The first integer
#' represents the chromosome number, and the second integer represents the
#' position on that chromosome. Th e encoded string can be used as a unique
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
#' @param chr Integer representing the chromosome number.
#' @param pos Integer representing the position on the chromosome.
#' @param max_chr Integer representing the maximum chromosome number expected.
#'   This is used to determine the fixed width for encoding the chromosome
#'   number. If left NULL, the maximum value of `chr` will be used.
#' @returns A base62 encoded string representing the combined chromosome and
#'   position.
#' @export
#' @examples
#' encoded_coords <- encode62(c(1,10,260),c(1,1000,1000000))
#' print(encoded_coords)
#' decoded_coords <- decode62(encoded_coords, max_chr = 260)
#' print(decoded_coords)

encode62 <- function(chr, pos, max_chr = NULL) {
  if (is.null(max_chr)) {
    max_chr <- max(chr)
  }
  # check number of digits needed for chr
  num_digits_chr <- ceiling(log(max_chr + 1, base = 62))
  # use the cpp encoding function
  encode_pair_vec_cpp(chr, pos, num_digits_chr)
}

#' Function to decode a base62 encoded string back into a pair of integers.
#' 
#' This function reverses the encoding performed by `encode_62`, taking a base62
#' encoded string and extracting the original chromosome number and position.
#' It uses the fixed width determined by the maximum chromosome number to
#' correctly separate the two components of the encoding.
#' @param encoded_str A base62 encoded string representing the combined
#'   chromosome and position
#' @param max_chr Integer representing the maximum chromosome number expected.
#'   This is used to determine the fixed width for decoding the chromosome
#'   number. If left NULL, an error will be raised as this information is needed
#'   for decoding.
#' @param num_digits_chr (Optional) Integer representing the number of characters
#'   used to encode the chromosome number. If provided, this will override the
#'   calculation based on `max_chr`.
#' @returns A data.frame with two columns: `chr` and `pos`, representing the
#'   decoded chromosome numbers and positions.
#' @export
#' @examples
#' encoded_coords <- encode62(c(1,10,260),c(1,1000,1000000))
#' print(encoded_coords)
#' decoded_coords <- decode62(encoded_coords, max_chr = 260)
#' print(decoded_coords)
decode62 <- function(encoded_str, max_chr = NULL, num_digits_chr = NULL) {
  # either max_chr or num_digits_chr must be provided
  if (is.null(max_chr) && is.null(num_digits_chr)) {
    stop("Either max_chr or num_digits_chr must be provided")
  }
  # check that they are not both defined
  if (!is.null(max_chr) && !is.null(num_digits_chr)) {
    stop("Both max_chr and num_digits_chr are provided; ",
    "only one can be given at a time")
  }
  if (!is.null(num_digits_chr)) {
    # use the provided num_digits_chr
    return(decode_pair_vec_cpp(encoded_str, num_digits_chr))
  }
  # check number of digits needed for chr
  num_digits_chr <- ceiling(log(max_chr + 1, base = 62))
  # use the cpp decoding function
  decode_pair_vec_cpp(encoded_str, num_digits_chr)
}
