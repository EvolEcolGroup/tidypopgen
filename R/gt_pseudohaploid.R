#' Set the ploidy of a `gen_tibble` to include pseudohaploids
#'
#' If a `gen_tibble` includes pseudohaploid data, its ploidy is set to -2 to
#' indicate that some individuals are coded as pseudohaploids. The ploidy of the
#' individuals is updated, with pseudohaploids set to 1 and diploids set to 2.
#' However, the dosages are not changed, meaning that pseudohaploids are still
#' coded as 0 or 2. If the `gen_tibble` is already set to pseudohaploid, running
#' gt_pseudohaploid will update the ploidy values again, if pseudohaploid
#' individuals have been removed then ploidy is reset to 2.
#'
#' @param x a `gen_tibble` object
#' @param test_n_loci the number of loci to test to determine if an individual
#' is pseudohaploid. If there are no heterozygotes in the first `test_n_loci`
#' loci, the individual is considered a pseudohaploid. If `NULL`, all loci are
#' tested.
#' @return a `gen_tibble` object with the ploidy set to -2 and the individual
#' ploidy values updated to 1 or 2.
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Detect pseudohaploids and set ploidy for the whole gen_tibble
#' example_gt <- example_gt %>% gt_pseudohaploid(test_n_loci = 3)
#'
#' # Ploidy is now set to -2
#' show_ploidy(example_gt)
#'
#' # Individual ploidy now varies between 1 (pseudohaploid) and 2 (diploid)
#' indiv_ploidy(example_gt)
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

  if (!"ploidy" %in% names(attr(x, "bigsnp")$fam)) {
    attr(x, "bigsnp")$fam$ploidy <- NA_integer_
  }

  attr(x$genotypes, "bigsnp")$fam$ploidy[.gt_bigsnp_rows(x)] <-
    identify_pseudohaploids(x, n_test = test_n_loci)

  if (min(attr(x$genotypes, "bigsnp")$fam$ploidy[.gt_bigsnp_rows(x)]) == 2) {
    # if all individuals are diploid, set ploidy to 2
    attr(x$genotypes, "ploidy") <- 2
  } else {
    # otherwise, set ploidy to -2
    attr(x$genotypes, "ploidy") <- -2
  }
  return(x)
}


#' Identify pseudohaploids
#'
#' Pseudohaploids are coded as all homozygotes; we find them by checking the
#' first 'n_test' loci and call a pseudohaploid if they have zero heterozygosity
#' (this is the same strategy employed in admixtools)
#' @param x the gen_tibble
#' @param n_test the number of loci being tested
#' @return a numeric vector of ploidy
#' @keywords internal
#' @noRd
identify_pseudohaploids <- function(x, n_test = 1000) {
  if (n_test > count_loci(x)) {
    n_test <- count_loci(x)
  }
  sub_x <- select_loci(x, .sel_arg = dplyr::all_of(c(1:n_test))) %>%
    dplyr::ungroup()
  het_obs <- indiv_het_obs(sub_x)
  ploidy <- rep(2, nrow(x))
  ploidy[het_obs == 0] <- 1
  return(ploidy)
}
