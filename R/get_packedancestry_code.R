#' Generate a packed ancestry code
#'
#' This function returns a matrix of raw values (which can be used in a file
#' backed matrix) that represent the packed ancestry code. The columns
#' are the packed bytes from 0 to 255, and the rows are the genotypes for the
#' four individuals represented by the byte. Each two bits give the number of
#' alleles, with the combination 11 representing a missing value.
#' @param minor A boolean indicating whether to use the minor allele or not.
#' *packedancestry* uses the major allele, but by default we convert to minor
#' alleles to be in line with PLINK and VCF.
#' @returns a matrix of raw values
#'
#' @keywords internal
#' @noRd
get_packedancestry_code <- function(minor = TRUE) { # TODO temp switch
  raw_vals <- as.raw(0:255)
  parsed_vals <- lapply(raw_vals, parse_2bit_groups)
  parsed_matrix <- do.call(rbind, parsed_vals)
  parsed_matrix <- t(parsed_matrix)
  if (minor) {
    # convert to minor alleles
    two_positions <- parsed_matrix == 2
    zero_positions <- parsed_matrix == 0
    parsed_matrix[two_positions] <- 0
    parsed_matrix[zero_positions] <- 2
  }
  storage.mode(parsed_matrix) <- "raw"

  return(parsed_matrix)
}

#' Parse a byte into 2-bit groups
#'
#' This function takes a byte (raw value) and extracts 4 groups of 2 bits each.
#' Following the convention by *packedancestry*, we start from the left (i.e the
#' higher value bits) and move to the right (i.e the lower value bits).
#' @param byte A raw value representing a byte (0-255).
#' @returns A vector of integers representing the 4 groups of 2 bits.
#' @keywords internal
#' @noRd
parse_2bit_groups <- function(byte) {
  if (!inherits(byte, "raw") || length(byte) != 1) {
    stop("Input must be a single raw value.")
  }

  byte_int <- as.integer(byte)

  # Extract 4 groups of 2 bits
  groups <- integer(4)
  for (i in 0:3) {
    mask <- bitwShiftL(3L, 2 * i)
    val <- bitwAnd(byte_int, mask)
    groups[i + 1] <- bitwShiftR(val, 2 * i)
  }

  # now reverse the order of the groups
  groups <- rev(groups)
  return(groups)
}
