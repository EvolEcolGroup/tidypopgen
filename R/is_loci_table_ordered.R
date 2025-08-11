#' Test if the loci table is ordered
#'
#' This functions checks that all SNPs in a chromosome are adjacent in the loci
#' table, and that positions are sorted within chromosomes.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param error_on_false logical, if `TRUE` an error is thrown if the loci are
#'   not ordered.
#' @param ignore_genetic_dist logical, if `TRUE` the physical position is not
#'   checked.
#' @param ... other arguments passed to specific methods.
#' @returns a logical vector defining which loci are transversions
#' @rdname is_loci_table_ordered
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% is_loci_table_ordered()
#'
is_loci_table_ordered <- function(
    .x,
    error_on_false = FALSE,
    ignore_genetic_dist = TRUE,
    ...) {
  UseMethod("is_loci_table_ordered", .x)
}

#' @export
#' @rdname is_loci_table_ordered
is_loci_table_ordered.tbl_df <- function(
    .x,
    error_on_false = FALSE,
    ignore_genetic_dist = TRUE,
    ...) {
  stopifnot_gen_tibble(.x)
  is_loci_table_ordered(
    .x$genotypes,
    error_on_false = error_on_false,
    ignore_genetic_dist = ignore_genetic_dist,
    ...
  )
}


#' @export
#' @rdname is_loci_table_ordered
is_loci_table_ordered.vctrs_bigSNP <- function(
    .x,
    error_on_false = FALSE,
    ignore_genetic_dist = TRUE,
    ...) {
  rlang::check_dots_empty()

  # check that within each chromosome positions are sorted
  if (
    any(unlist(
      show_loci(.x) %>%
        group_by(.data$chr_int) %>% # nolint
        group_map(~ is.unsorted(.x$position))
    ))
  ) {
    if (error_on_false) {
      stop("Your loci are not sorted within chromosomes")
    } else {
      return(FALSE)
    }
  }

  # check that within each chromosome positions are unique
  if (isFALSE(find_duplicated_loci(
    .x,
    error_on_false = error_on_false,
    list_duplicates = FALSE
  ))) {
    return(FALSE)
  }

  # check that all positions in a chromosome are adjacent
  if (any(duplicated(rle(show_loci(.x)$chr_int)$values))) {
    if (error_on_false) {
      stop("All SNPs in a chromosome should be adjacent in the loci table")
    } else {
      return(FALSE)
    }
  }

  # check genetic distance
  if (!ignore_genetic_dist) {
    if (
      any(unlist(
        show_loci(.x) %>%
          group_by(.data$chr_int) %>% # nolint
          group_map(~ is.unsorted(.x$genetic_dist))
      ))
    ) {
      if (error_on_false) {
        stop("Your genetic distances are not sorted within chromosomes")
      } else {
        return(FALSE)
      }
    }
    if (
      any(unlist(
        show_loci(.x) %>%
          group_by(.data$chr_int) %>% # nolint
          group_map(~ duplicated(.x$genetic_dist))
      ))
    ) {
      if (all(show_loci(.x)$genetic_dist == 0) && error_on_false) {
        message("Your genetic distances have been set to 0")
        return(TRUE)
      } else if (all(show_loci(.x)$genetic_dist == 0) && !error_on_false) {
        return(TRUE)
      } else if (!all(show_loci(.x)$genetic_dist == 0) && error_on_false) {
        stop("Your loci table contains duplicated genetic distances")
      } else {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}


#' Find duplicates in the loci table
#'
#' This function finds duplicated SNPs by checking the positions within each
#' chromosome. It can return a list of duplicated SNPs or a logical value
#' indicating whether there are any duplicated loci.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param error_on_false logical, if `TRUE` an error is thrown if duplicated
#'   loci are found.
#' @param list_duplicates logical, if `TRUE` returns duplicated SNP names.
#' @param ... other arguments passed to specific methods.
#' @returns if `list_duplicates` is TRUE, returns a vector of duplicated loci.
#'   If `list_duplicates` is FALSE, returns a logical value indicating whether
#'   there are any duplicated loci. If `error_on_false` is TRUE and there are
#'   duplicates, an error is thrown.
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' show_loci(example_gt) <- test_loci <- data.frame(
#'   big_index = c(1:6),
#'   name = paste0("rs", 1:6),
#'   chromosome = paste0("chr", c(1, 1, 1, 1, 1, 1)),
#'   position = as.integer(c(3, 3, 5, 65, 343, 46)),
#'   genetic_dist = as.double(rep(0, 6)),
#'   allele_ref = c("A", "T", "C", "G", "C", "T"),
#'   allele_alt = c("T", "C", NA, "C", "G", "A"),
#'   chr_int = rep(1, 6)
#' )
#'
#' show_loci(example_gt)
#'
#' # Find which loci are duplicated
#' example_gt %>% find_duplicated_loci()
find_duplicated_loci <- function(.x,
                                 error_on_false = FALSE,
                                 list_duplicates = TRUE,
                                 ...) {
  if (
    any(unlist(
      show_loci(.x) %>%
        group_by(.data$chr_int) %>% # nolint
        group_map(~ duplicated(.x$position))
    ))
  ) {
    if (list_duplicates) {
      show_loci(.x) %>%
        group_by(.data$chr_int) %>% # nolint
        group_map(~ .x[duplicated(.x$position) | duplicated(.x$position, fromLast = TRUE), ]$name) %>% # nolint
        unlist(use.names = FALSE)
    } else if (error_on_false) {
      stop("Your loci table contains duplicates")
    } else {
      return(FALSE)
    }
  } else {
    return(TRUE)
  }
}
