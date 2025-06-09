#' Flip the strand of alleles
#'
#' This function takes a vector of bases and recodes them.
#'
#' Note that it does not check that all bases are valid. Also, missing alleles
#' "m" are kept as "m".
#'
#' @param bases a vector of characters. The only valid characters are "A", "T",
#' "C" and "G".
#' @returns a vector of recoded bases
#' @keywords internal
#' @noRd
flip <- function(bases) {
  to_a <- bases == "T"
  to_t <- bases == "A"
  to_c <- bases == "G"
  to_g <- bases == "C"
  bases[to_a] <- "A"
  bases[to_t] <- "T"
  bases[to_c] <- "C"
  bases[to_g] <- "G"
  return(bases) # nolint
}
