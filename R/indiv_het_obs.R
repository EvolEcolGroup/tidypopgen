#' Estimate individual observed heterozygosity
#'
#' Estimate observed heterozygosity (H_obs) for each individual (i.e. the
#' frequency of loci that are heterozygous in an individual). A locus is
#' heterozygous whenever the dosage of the alternate allele is neither zero
#' nor the individual's own ploidy, so this works for diploid, polyploid, and
#' mixed-ploidy [`gen_tibble`] objects alike (each individual is evaluated
#' against its own ploidy, so individuals of differing ploidy can be mixed
#' freely).
#'
#' On a grouped [`gen_tibble`], heterozygosity is computed for all the
#' individuals in a single pass over the genotypes (it does not depend on the
#' groups), and returned in the order of the rows. Within a grouped
#' `mutate()` or `summarise()`, `dplyr` calls `indiv_het_obs()` on the
#' genotypes of each group in turn: this gives the same values, but reads the
#' genotypes once for each group. With many groups and large datasets,
#' compute it on the whole grouped tibble instead (e.g.
#' `grouped_gt$het_obs <- indiv_het_obs(grouped_gt)`), or before grouping.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`] (grouped or not).
#' @param as_counts logical, if `TRUE`, return a matrix with two columns:
#'  the number of heterozygotes and the number of missing values for each
#'  individual. These quantities can be useful to compute more complex
#'  quantities.
#' @param ... currently unused.
#' @returns either:
#'  - a vector of heterozygosities, one per individuals in the [`gen_tibble`]
#'  - a matrix with two columns, where the first is the number of heterozygous
#'    loci for each individual and the second is the number of missing values
#'    for each individual
#' @rdname indiv_het_obs
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% indiv_het_obs()
#'
#' # For observed heterozygosity as counts:
#' example_gt %>% indiv_het_obs(as_counts = TRUE)
#'
indiv_het_obs <- function(.x, as_counts = FALSE, ...) {
  UseMethod("indiv_het_obs", .x)
}

#' @export
#' @rdname indiv_het_obs
indiv_het_obs.tbl_df <- function(.x, as_counts = FALSE, ...) {
  stopifnot_has_genotypes(.x)
  # extract the column and hand it over to its method
  indiv_het_obs(.x$genotypes, as_counts = as_counts, ...)
}

#' @export
#' @rdname indiv_het_obs
indiv_het_obs.grouped_gen_tbl <- function(.x, as_counts = FALSE, ...) {
  stopifnot_has_genotypes(.x)
  # heterozygosity does not depend on the groups: compute it for all the
  # individuals in a single pass over the genotypes (the whole column), in
  # the order of the rows
  indiv_het_obs(.x$genotypes, as_counts = as_counts, ...)
}

# a grouped data frame that has lost the class of its tibble (e.g. a grouped
# gen_tibble returned by summarise()): restore the class and dispatch again
#' @export
#' @noRd
indiv_het_obs.grouped_df <- function(.x, ...) {
  indiv_het_obs(restore_genotypes_tbl(.x, "indiv_het_obs"), ...)
}

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
indiv_het_obs.data.frame <- function(.x, ...) {
  indiv_het_obs(restore_genotypes_tbl(.x, "indiv_het_obs"), ...)
}

#' @export
#' @rdname indiv_het_obs
indiv_het_obs.vctrs_bigSNP <- function(.x, as_counts = FALSE, ...) {
  rlang::check_dots_empty()
  stopifnot_no_pseudohaploid(.x)
  # get the FBM
  X <- attr(.x, "fbm") # nolint
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)

  # use a fast (and simpler) diploid/pseudohaploid version
  # (no per-individual lookup or comparison against ploidy), and only pay
  # for the more general, ploidy-aware kernel when the data actually need
  # it, so the common diploid case is not slowed down.
  if (is_diploid_only(.x) || is_pseudohaploid(.x)) {
    # returns a matrix of 2 rows (count_1,count_na) and n_individuals columns
    count_1_na <- function(BM, ind, rows_to_keep) { # nolint
      gt_ind_hetero(
        BM = BM,
        rowInd = rows_to_keep,
        colInd = ind,
        ncores = 1 # n_cores, I have not seen any improvement with n_cores > 1
      )
    }
    this_col_1_na <- bigstatsr::big_apply(
      X,
      a.FUN = count_1_na,
      ind = attr(.x, "loci")$big_index,
      a.combine = "plus",
      rows_to_keep = rows_to_keep
    )
  } else {
    # ploidy aware version
    # heterozygosity is an individual-level property, so ploidy is looked
    # up per individual; this works even for mixed-ploidy gen_tibbles
    ploidy <- indiv_ploidy(.x)
    count_1_na_poly <- function(BM, ind, rows_to_keep, ploidy) { # nolint
      gt_ind_hetero_polyploid(
        BM = BM,
        rowInd = rows_to_keep,
        colInd = ind,
        ploidy = ploidy,
        ncores = 1 # n_cores, I have not seen any improvement with n_cores > 1
      )
    }
    this_col_1_na <- bigstatsr::big_apply(
      X,
      a.FUN = count_1_na_poly,
      ind = attr(.x, "loci")$big_index,
      a.combine = "plus",
      rows_to_keep = rows_to_keep,
      ploidy = ploidy
    )
  }
  if (!as_counts) {
    # use the number of selected loci, not the number of columns in the FBM
    n_loci <- length(attr(.x, "loci")$big_index)
    return(this_col_1_na[1, ] / (n_loci - this_col_1_na[2, ]))
  } else {
    this_col_1_na <- t(this_col_1_na)
    colnames(this_col_1_na) <- c("het_n", "na_n")
    return(this_col_1_na)
  }
}
