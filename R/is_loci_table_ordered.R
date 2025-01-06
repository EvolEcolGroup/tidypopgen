#' Test if the loci table is ordered
#'
#' This functions checks that all SNPs in a chromosome are adjacent in the loci table,
#' and that positions are sorted within chromosomes.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param error_on_false logical, if `TRUE` an error is thrown if the loci
#' are not ordered.
#' @param ignore_genetic_dist logical, if `TRUE` the physical position is not checked.
#' @param ... other arguments passed to specific methods.
#' @returns a logical vector defining which loci are transversions
#' @rdname is_loci_table_ordered
#' @export
is_loci_table_ordered <- function(.x, error_on_false = FALSE, ignore_genetic_dist = TRUE, ...) {
  UseMethod("is_loci_table_ordered", .x)
}

#' @export
#' @rdname is_loci_table_ordered
is_loci_table_ordered.tbl_df <- function(.x, error_on_false = FALSE, ignore_genetic_dist = TRUE, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  is_loci_table_ordered(.x$genotypes, error_on_false = error_on_false,
                        ignore_genetic_dist = ignore_genetic_dist, ...)
}


#' @export
#' @rdname is_loci_table_ordered
is_loci_table_ordered.vctrs_bigSNP <- function(.x, error_on_false = FALSE, ignore_genetic_dist = TRUE, ...) {
  rlang::check_dots_empty()

  # check that within each chromosome positions are sorted
  if (any(unlist(show_loci(.x) %>%
                 group_by(.data$chr_int) %>%
                 group_map(~ is.unsorted(.x$position))))){
    if (error_on_false){
      stop("Your loci are not sorted within chromosomes")
    } else {
      return(FALSE)
    }
  }

  # check that within each chromosome positions are unique
  if (any(unlist(show_loci(.x) %>%
                 group_by(.data$chr_int) %>%
                 group_map(~ duplicated(.x$position))))){
    if (error_on_false){
      stop("Your loci are not sorted within chromosomes")
    } else {
      return(FALSE)
    }
  }

  # check that all positions in a chromosome are adjacent
  if(any(duplicated(rle(show_loci(.x)$chr_int)$values))){
    if (error_on_false){
      stop("All SNPs in a chromosome should be adjacent in the loci table")
    } else {
      return(FALSE)
    }
  }

  # check genetic distance
  if (!ignore_genetic_dist){
    if (any(unlist(show_loci(.x) %>%
                   group_by(.data$chr_int) %>%
                   group_map(~ is.unsorted(.x$genetic_dist))))){
      if (error_on_false){
        stop("Your genetic distances are not sorted within chromosomes")
      } else {
        return(FALSE)
      }
    }
  }

  return(TRUE)
}

