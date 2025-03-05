#' Estimate individual observed heterozygosity
#'
#' Estimate observed heterozygosity (H_obs) for each individual (i.e. the
#' frequency of loci that are heterozygous in an individual).
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a vector of heterozygosities, one per individuals in the
#'   [`gen_tibble`]
#' @rdname indiv_het_obs
#' @export
indiv_het_obs <- function(.x, ...) {
  UseMethod("indiv_het_obs", .x)
}

#' @export
#' @rdname indiv_het_obs
indiv_het_obs.tbl_df <- function(.x, ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  indiv_het_obs(.x$genotypes, ...)
}

#' @export
#' @rdname indiv_het_obs
indiv_het_obs.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  # get the FBM
  X <- attr(.x, "bigsnp")$genotypes # nolint
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # col means for submatrix (all rows, only some columns)

  # returns a matrix of 2 rows (count_1,count_na) and n_individuals columns
  col_1_na <- function(X, ind, rows_to_keep) { # nolint
    count_1 <- function(a) {
      sum(a == 1, na.rm = TRUE)
    }
    res <- apply(X[rows_to_keep, ind], 1, count_1)
    count_na <- function(a) {
      sum(is.na(a))
    }
    res <- rbind(res, apply(X[rows_to_keep, ind], 1, count_na))
  }

  # count heterozygotes and nas in one go
  this_col_1_na <- bigstatsr::big_apply(X,
    a.FUN = col_1_na,
    ind = attr(.x, "loci")$big_index,
    a.combine = "plus", rows_to_keep = rows_to_keep
  )
  this_col_1_na[1, ] / (ncol(X) - this_col_1_na[2, ])
}

# #' @export #nolint start
# #' @rdname indiv_het_obs
# indiv_het_obs.grouped_df <- function(.x, ...){
#   .x %>% mutate(indiv_het_obs = indiv_het_obs(.data$genotypes)) %>% summarise(het_obs = mean(indiv_het_obs))
# } #nolint end
