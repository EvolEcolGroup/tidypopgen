#' Set the ploidy of a `gen_tibble` to include pseudohaploids
#' 
#' The ploidy of a `gen_tibble` is set to -2, to indicate that some individuals
#' are coded as pseudohaploids. The ploidy of the individuals is updated, with
#' pseudohaploids set to 1 and diploids set to 2. However, the dosages are not
#' changed, meaning that pseudohaploids are still coded as 0 or 2.
#' 
#' @param x a `gen_tibble` object
#' @param test_n_loci the number of loci to test to determine if an individual
#' is pseudohaploid. If there are no heterozygotes in the first `test_n_loci`
#' loci, the individual is considered a pseudohaploid. If `NULL`, all loci are
#' tested.
#' @return a `gen_tibble` object with the ploidy set to -2 and the individual
#' ploidy values updated to 1 or 2.
#' @export

gt_pseudohaploid <- function(x, test_n_loci = 10000) {
  # check that the input is a gen_tibble
  stopifnot_gen_tibble(x)
  
  # if this is already set to pseudohaploid, reset it to diploid to reprocess it
  if (attr(x$genotypes, "ploidy") == -2) {
    attr(x$genotypes, "ploidy") <- 2
  }

  # get the number of loci
  n_loci <- count_loci(x)
  
  # if test_n_loci is NULL, set it to n_loci
  if (is.null(test_n_loci)) {
    test_n_loci <- n_loci
  }
  # if test_n_loci greater than n_loci, set it to n_loci
  test_n_loci <- min(test_n_loci, n_loci)
  
  if (!"ploidy" %in% names(attr(x, "bigsnp")$fam)){
    attr(x, "bigsnp")$fam$ploidy <- NA_integer_
  }

  attr(x$genotypes, "bigsnp")$fam$ploidy[.gt_bigsnp_rows(x)] <-
    identify_pseudohaploids (x, n_test = test_n_loci)
  
  # set ploidy to -2
  attr(x$genotypes, "ploidy") <- -2
  
  
  return(x)
}
