#' Estimate missingness at each locus
#'
#' Estimate the rate of missingness at each locus.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param as_counts boolean defining whether the count of NAs (rather than the rate)
#' should be returned. It defaults to FALSE (i.e. rates are returned by default).
#' @param n_cores number of cores to be used, it defaults to [bigstatsr::nb_cores()]
#' @param ... other arguments passed to specific methods.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_missingness
#' @export
loci_missingness <- function(.x, as_counts = FALSE, ...) {
  UseMethod("loci_missingness", .x)
}

#' @export
#' @rdname loci_missingness
loci_missingness.tbl_df <- function(.x, as_counts = FALSE, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_missingness(.x$genotypes, as_counts = as_counts, ...)
}


#' @export
#' @rdname loci_missingness
loci_missingness.vctrs_bigSNP <- function(.x, as_counts = FALSE, ...) {
  rlang::check_dots_empty()
#  stopifnot_diploid(.x)
  # get the FBM
  geno_fbm <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep)>1){
    n_na <- bigstatsr::big_counts(geno_fbm, ind.row = rows_to_keep,
                          ind.col = attr(.x,"loci")$big_index)
    n_na <- n_na[nrow(n_na),] # this should work also with polyploids
    if (!as_counts){
      n_na <- n_na/length(rows_to_keep)
    }
  } else { # if we have a single individual
    n_na <-geno_fbm[rows_to_keep,attr(.x,"loci")$big_index]
  }
  n_na
}

#' @export
#' @rdname loci_missingness
loci_missingness.grouped_df <- function(.x, as_counts = FALSE,
                                        n_cores = bigstatsr::nb_cores(), ...) {
  rlang::check_dots_empty()
    geno_fbm <- .gt_get_bigsnp(.x)$genotypes

    na_mat <- gt_grouped_missingness(BM = geno_fbm,rowInd = .gt_bigsnp_rows(.x),
                                            colInd = .gt_bigsnp_cols(.x),
                                            groupIds = dplyr::group_indices(.x)-1,
                                            ngroups = max(dplyr::group_indices(.x)),
                                            ncores = n_cores)
    group_sizes <- tally(.x) %>% dplyr::pull(dplyr::all_of("n"))
    if (!as_counts){
     na_mat <- sweep(na_mat, MARGIN=2, STATS=group_sizes,FUN="/")
    }

    # return a list to mimic a group_map
    lapply(seq_len(ncol(na_mat)), function(i) na_mat[,i])
}

