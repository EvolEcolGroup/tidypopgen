#' Estimate allele frequencies at each each locus
#'
#' Allele frequencies can be estimates as minimum allele frequencies (MAF) with
#' `loci_maf()` or the frequency of the alternate allele (with `loci_alt_freq()`).
#' The latter are in line with the genotypes matrix (e.g. as extracted by
#'  [`show_loci()`]). Most users will be in interested in the MAF, but the
#'  raw frequencies might be useful when computing aggregated statistics.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param n_cores number of cores to be used, it defaults to [bigstatsr::nb_cores()]
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_alt_freq
#' @export
loci_alt_freq <- function(.x, ...) {
  UseMethod("loci_alt_freq", .x)
}

#' @export
#' @rdname loci_alt_freq
loci_alt_freq.tbl_df <- function(.x, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_alt_freq(.x$genotypes)
}


#' @export
#' @rdname loci_alt_freq
loci_alt_freq.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  #stopifnot_diploid(.x)
  # if we have diploid
  if (is_diploid_only(.x)){
    loci_alt_freq_diploid(.x)
  } else {
    loci_alt_freq_polyploid(.x)
  }
}

#' @export
#' @rdname loci_alt_freq
loci_alt_freq.grouped_df <- function(.x, n_cores = bigstatsr::nb_cores(), ...) {
  rlang::check_dots_empty()
  if (is_diploid_only(.x)){
    geno_fbm <- .gt_get_bigsnp(.x)$genotypes

    freq_mat <- gt_grouped_alt_freq_diploid(BM = geno_fbm,rowInd = .gt_bigsnp_rows(.x),
                  colInd = .gt_bigsnp_cols(.x),
                  groupIds = dplyr::group_indices(.x)-1,
                  ngroups = max(dplyr::group_indices(.x)),
                  ncores = n_cores)$freq_alt
    # return a list to mimic a group_map
    lapply(seq_len(ncol(freq_mat)), function(i) freq_mat[,i])
  } else {
    # TODO this is seriously inefficient, we should replace it with a cpp function
    group_map(.x, .f=~loci_alt_freq(.x,, ...))
  }


}

#' @rdname loci_alt_freq
#' @export
loci_maf <- function(.x, ...) {
  UseMethod("loci_maf", .x)
}

#' @export
#' @rdname loci_alt_freq
loci_maf.tbl_df <- function(.x, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_maf(.x$genotypes, ...)
}

#' @export
#' @rdname loci_alt_freq
loci_maf.vctrs_bigSNP <- function(.x, ...) {
  freq <- loci_alt_freq(.x,)
  freq[freq>0.5 & !is.na(freq)] <- 1 - freq[freq>0.5 & !is.na(freq)]
  freq
}

#' @export
#' @rdname loci_alt_freq
loci_maf.grouped_df <- function(.x, n_cores = bigstatsr::nb_cores(), ...) {
  rlang::check_dots_empty()
  if (is_diploid_only(.x)){
    geno_fbm <- .gt_get_bigsnp(.x)$genotypes

    freq_mat <- gt_grouped_alt_freq_diploid(BM = geno_fbm,rowInd = .gt_bigsnp_rows(.x),
                                            colInd = .gt_bigsnp_cols(.x),
                                            groupIds = dplyr::group_indices(.x)-1,
                                            ngroups = max(dplyr::group_indices(.x)),
                                            ncores = n_cores)$freq_alt
    freq_mat[freq_mat>0.5 & !is.na(freq_mat)] <- 1 - freq_mat[freq_mat>0.5 & !is.na(freq_mat)]
    # return a list to mimic a group_map
    lapply(seq_len(ncol(freq_mat)), function(i) freq_mat[,i])
  } else { # the polyploid case
    # TODO this is seriously inefficient, we should replace it with a cpp function
    group_map(.x, .f=~loci_maf(.x,, ...))
  }
}

# function to estimate frequencies for diploid
loci_alt_freq_diploid <- function(.x){
  # get the FBM
  geno_fbm <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep)>1){
    # create function to use in big_apply
    big_sub_counts <- function (X, ind, rows_to_keep) {
      col_counts <- bigstatsr::big_counts(X,
                                          ind.row = rows_to_keep,
                                          ind.col = ind)
      means_from_counts <- function(x){
        (x[2]+x[3]*2)/((x[1]+x[2]+x[3])*2)
      }
      freq_sub <- apply(col_counts, 2, means_from_counts)
      freq_sub
    }
    freq <- bigstatsr::big_apply(geno_fbm, a.FUN = big_sub_counts,
                                 rows_to_keep = rows_to_keep,
                                 ind=attr(.x,"loci")$big_index,
                                 ncores = 1, # we only use 1 cpu, we let openMP use multiple cores
                                 block.size = bigstatsr::block_size(attr(.x,"loci")$big_index, 1),
                                 a.combine = 'c')

    # col_counts <- bigstatsr::big_counts(geno_fbm,
    #                                     ind.row = rows_to_keep,
    #                                     ind.col = attr(.x,"loci")$big_index)
    # means_from_counts <- function(x){
    #   (x[2]+x[3]*2)/((x[1]+x[2]+x[3])*2)
    # }
    # freq <- apply(col_counts, 2, means_from_counts)
  } else { # if we have a single individual
    freq <-geno_fbm[rows_to_keep,attr(.x,"loci")$big_index] /2
  }
  freq
}

loci_alt_freq_polyploid <- function(.x, ...){
  warning("this function still needs a proper unit test!!! It assumes alleles are the unit of observation")
  # get the FBM
  geno_fbm <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  ploidy_by_indiv <- indiv_ploidy(.x)
  if (length(rows_to_keep)>1){
    # col means for submatrix (all rows, only some columns)
    col_sums_na <- function(X, ind, rows_to_keep, ploidy_by_indiv) {

      res <- colSums(X[rows_to_keep, ind], na.rm=TRUE)
      col_na <- function(a, ploidy_by_indiv){sum(is.na(a)*ploidy_by_indiv)}
      res <- cbind(res,apply(X[rows_to_keep,ind],2,col_na, ploidy_by_indiv = ploidy_by_indiv))
      res
    }
    col_sums_na_mat <- bigstatsr::big_apply(geno_fbm, a.FUN = col_sums_na,
                                 rows_to_keep = rows_to_keep,
                                 ind=attr(.x,"loci")$big_index,
                                 ploidy_by_indiv = ploidy_by_indiv,
                                 a.combine = 'rbind')
    # now get frequency accounting for missing values
    col_sums_na_mat[,1]/(sum(ploidy_by_indiv) - col_sums_na_mat[,2])
  } else { # if we have a single individual
    geno_fbm[rows_to_keep,attr(.x,"loci")$big_index]/ploidy_by_indiv
  }
}
