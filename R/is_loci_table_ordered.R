#' Test if the loci table is ordered
#'
#' This functions checks that all SNPs in a chromosome are adjacent in the loci
#' table, and that positions (and, unless `ignore_genetic_dist = TRUE`,
#' genetic distances) are sorted within chromosomes. Missing positions (or
#' missing genetic distances, if they are checked) are an error, as the order
#' of those loci is unknown.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param error_on_false logical, if `TRUE` an error is thrown if the loci are
#'   not ordered.
#' @param ignore_genetic_dist logical, if `TRUE` (the default) only the
#'   positions are checked, and genetic distances are ignored (they can then
#'   be missing).
#' @param ... other arguments passed to specific methods.
#' @returns `TRUE` if the loci are ordered, `FALSE` otherwise
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
  ...
) {
  UseMethod("is_loci_table_ordered", .x)
}

#' @export
#' @rdname is_loci_table_ordered
is_loci_table_ordered.tbl_df <- function(
  .x,
  error_on_false = FALSE,
  ignore_genetic_dist = TRUE,
  ...
) {
  stopifnot_has_genotypes(.x)
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
  ...
) {
  rlang::check_dots_empty()

  loci <- show_loci(.x)

  # loci with a missing position (or genetic distance, if it is checked) can
  # not be ordered
  if (anyNA(loci$position)) {
    stop(
      sum(is.na(loci$position)), " loci have a missing position, so the ",
      "order of the loci table can not be checked"
    )
  }
  if (!ignore_genetic_dist && anyNA(loci$genetic_dist)) {
    stop(
      sum(is.na(loci$genetic_dist)), " loci have a missing genetic distance; ",
      "use ignore_genetic_dist = TRUE to check the order of positions only"
    )
  }

  # check that within each chromosome positions are sorted (each check is a
  # single pass in C++, which does not copy or split the columns of the loci
  # table)
  position_order <- any_unsorted_by_chromosome(loci$position, loci$chromosome)
  if (position_order[["unsorted"]]) {
    if (error_on_false) {
      stop("Your loci are not sorted within chromosomes")
    } else {
      return(FALSE)
    }
  }

  # check that within each chromosome positions are unique; as positions are
  # now known to be sorted, duplicated positions are equal consecutive
  # positions within a chromosome
  if (position_order[["ties"]]) {
    if (error_on_false) {
      stop("Your loci table contains duplicates")
    } else {
      return(FALSE)
    }
  }

  # check that all positions in a chromosome are adjacent, i.e. that no
  # chromosome has more than one run of loci (as with rle(), each missing
  # chromosome is a run of its own, so more than one missing chromosome fails)
  if (chromosomes_not_adjacent(loci$chromosome)) {
    if (error_on_false) {
      stop("All SNPs in a chromosome should be adjacent in the loci table")
    } else {
      return(FALSE)
    }
  }

  # check genetic distance (if they are sorted across the whole table, e.g.
  # all zero, they are sorted within each chromosome, and we skip grouping
  # the loci by chromosome; missing values were ruled out above)
  if (!ignore_genetic_dist) {
    if (is.unsorted(loci$genetic_dist) &&
      any_unsorted_by_chromosome(
        loci$genetic_dist,
        loci$chromosome
      )[["unsorted"]]) {
      if (error_on_false) {
        stop("Your genetic distances are not sorted within chromosomes")
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
#' @returns If `list_duplicates` is TRUE, returns a character vector of
#'   duplicated loci names (character(0) when none). If `list_duplicates` is
#'   FALSE, returns TRUE when no duplicates exist and FALSE when duplicates are
#'   present. If `error_on_false` is TRUE and duplicates exist, an error is
#'   thrown.
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
#'   allele_alt = c("T", "C", NA, "C", "G", "A")
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
  loci <- show_loci(.x)
  chr_int <- cast_chromosome_to_int(loci$chromosome)
  # pairs of chromosome and position (this avoids splitting the whole loci
  # table by chromosome)
  chr_pos <- vctrs::data_frame(chromosome = chr_int, position = loci$position)
  if (vctrs::vec_duplicate_any(chr_pos)) {
    if (list_duplicates) {
      # flag all loci that share chromosome and position with another locus
      is_dup <- vctrs::vec_duplicate_detect(chr_pos)
      # return their names ordered by chromosome (as when grouping by
      # chromosome), and then by their order in the loci table (the radix
      # sort is stable, and puts missing chromosomes last as group_by does)
      dup_index <- which(is_dup)
      dup_index <- dup_index[order(chr_int[dup_index], method = "radix")]
      loci$name[dup_index]
    } else if (error_on_false) {
      stop("Your loci table contains duplicates")
    } else {
      return(FALSE)
    }
  } else {
    if (list_duplicates) {
      return(character(0))
    } else {
      return(TRUE)
    }
  }
}


#' Test whether values are unsorted within any chromosome
#'
#' Equivalent to splitting `x` by `cast_chromosome_to_int()` of `chromosome`
#' and checking whether any of the groups is unsorted (as [is.unsorted()], so a
#' group of more than one value with missing values gives NA), but in a single
#' pass in C++, without copying or splitting `x`. Chromosomes that cast to the
#' same integer (e.g. "chr1" and "1") are a single group, as when sorting the
#' loci table in [gt_order_loci()].
#' @param x a numeric vector (e.g. the positions of the loci)
#' @param chromosome the chromosomes of the loci (a factor, as in the loci
#'   table)
#' @returns a named logical vector: `unsorted` is TRUE, FALSE or NA as
#'   described above, and `ties` is TRUE if two consecutive non-missing values
#'   of a chromosome are equal (when no chromosome is unsorted, this means that
#'   some pairs of chromosome and value are duplicated)
#' @keywords internal
#' @noRd
any_unsorted_by_chromosome <- function(x, chromosome) {
  # integer chromosomes (chromosomes_as_int) are already cast to integers,
  # each with its own group
  if (is.integer(chromosome)) {
    groups <- chromosome_group_codes(chromosome)
    return(unsorted_by_chromosome_cpp(x, groups$codes, seq_len(groups$n)))
  }
  if (!is.factor(chromosome)) {
    chromosome <- as.factor(chromosome)
  }
  # cast the levels rather than each locus; levels with the same integer
  # chromosome get the same group (the index of the first of those levels),
  # and levels that cast to NA are grouped with the missing chromosomes
  levels_int <- cast_chromosome_to_int(levels(chromosome))
  level_group <- match(levels_int, levels_int, incomparables = NA)
  unsorted_by_chromosome_cpp(x, chromosome, level_group)
}


#' Test whether any chromosome is split into several runs of loci
#'
#' Equivalent to finding the runs of chromosomes with [rle()] and checking
#' whether any chromosome has more than one run (so, as `rle()` does not
#' merge missing values, more than one missing chromosome gives TRUE), but in
#' a single pass in C++.
#' @param chromosome the chromosomes of the loci (a factor, as in the loci
#'   table)
#' @returns TRUE or FALSE
#' @keywords internal
#' @noRd
chromosomes_not_adjacent <- function(chromosome) {
  groups <- chromosome_group_codes(chromosome)
  chromosomes_not_adjacent_cpp(groups$codes, groups$n)
}

# The chromosomes of the loci as group codes 1..n (NA for missing
# chromosomes), for the C++ helpers that take the codes of a factor: the codes
# of a factor (with n its number of levels), or, for integer chromosomes
# (chromosomes_as_int), the integers themselves if they are positive and not
# much larger than the number of loci (so that the helpers' tables of n
# groups stay small), or else their index among the unique values (which are
# then returned as `values`, to map the groups back to the integers). Other
# vectors are turned into a factor.
chromosome_group_codes <- function(chromosome) {
  if (is.integer(chromosome) && !is.factor(chromosome)) {
    if (length(chromosome) == 0 || all(is.na(chromosome))) {
      return(list(codes = chromosome, n = 0L))
    }
    code_range <- range(chromosome, na.rm = TRUE)
    if (code_range[1] >= 1 &&
          code_range[2] <= max(2 * length(chromosome), 1e6)) {
      return(list(codes = chromosome, n = code_range[2]))
    }
    values <- unique(chromosome)
    return(list(
      codes = match(chromosome, values, incomparables = NA),
      n = length(values),
      values = values
    ))
  }
  if (!is.factor(chromosome)) {
    chromosome <- as.factor(chromosome)
  }
  list(codes = chromosome, n = nlevels(chromosome))
}
