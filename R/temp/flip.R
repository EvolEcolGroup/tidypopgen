#' Flip the strand of alleles
#'
#' This function takes a vector of bases and recodes them.
#'
#' Note that it does not check that all bases are valid. Also, missing alleles
#' "m" are kept as "m".
#'
#' @param bases a vector of characters. The only valid characters are "a", "t",
#' "c" and "g".
#' @returns a vector of recoded bases
#' @keywords internal

flip <-function(bases){
  to_a <- bases=="t"
  to_t <- bases=="a"
  to_c <- bases=="g"
  to_g <- bases=="c"
  bases[to_a] <- "a"
  bases[to_t] <- "t"
  bases[to_c] <- "c"
  bases[to_g] <- "g"
  return(bases)
}
