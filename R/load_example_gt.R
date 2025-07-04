#' Load example gen_tibble
#'
#' This function creates a `gen_tibble` object for use in examples in
#' documentation.
#' @param type a character string indicating the type of `gen_tibble` to create:
#'   - "gen_tbl": a basic gen_tibble with genotype data and metadata
#'   - "grouped_gen_tbl": same as "gen_tbl" but grouped by population
#'   - "grouped_gen_tbl_sf": adds spatial features (longitude/latitude)
#'      and groups by population
#'   - "gen_tbl_sf": adds spatial features without grouping
#' @returns an example object of the class `gen_tbl`.
#' @rdname load_example_gt
#' @export
#' @examples
#' # This function creates an example gen_tibble object
#' example_gt <- load_example_gt("gen_tbl")
load_example_gt <- function(type = c(
                              "gen_tbl", # nolint
                              "grouped_gen_tbl",
                              "grouped_gen_tbl_sf",
                              "gen_tbl_sf"
                            )) {
  type <- match.arg(type)
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, NA, 0, 0),
    c(2, NA, 0, 0, 1, 1),
    c(1, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(0, 0, 0, 0, NA, 1),
    c(0, 1, 1, 0, 1, NA)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e", "f", "g"),
    population = c("pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3")
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(3, 5, 65, 343, 23, 456)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  # Add spatial coordinates if needed
  if (type %in% c("grouped_gen_tbl_sf", "gen_tbl_sf")) {
    test_indiv_meta$longitude <- c(0, 0, 2, 2, 0, 2, 2)
    test_indiv_meta$latitude <- c(51, 51, 49, 49, 51, 41, 41)
  }

  # Create the base gen_tibble
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  # Add spatial features if needed
  if (type %in% c("grouped_gen_tbl_sf", "gen_tbl_sf")) {
    test_gt <- gt_add_sf(
      x = test_gt,
      coords = c("longitude", "latitude")
    )
  }

  # Add grouping if needed
  if (type %in% c("grouped_gen_tbl", "grouped_gen_tbl_sf")) {
    test_gt <- test_gt %>% group_by(.data$population)
  }
  return(test_gt)
}
