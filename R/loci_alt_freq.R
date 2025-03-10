#' Estimate allele frequencies at each each locus
#'
#' Allele frequencies can be estimates as minimum allele frequencies (MAF) with
#' `loci_maf()` or the frequency of the alternate allele (with
#' `loci_alt_freq()`). The latter are in line with the genotypes matrix (e.g. as
#' extracted by [`show_loci()`]). Most users will be in interested in the MAF,
#' but the raw frequencies might be useful when computing aggregated statistics.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param block_size maximum number of loci read at once.
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_alt_freq
#' @export
loci_alt_freq <- function(.x, n_cores, block_size, ...) {
  UseMethod("loci_alt_freq", .x)
}

#' @export
#' @rdname loci_alt_freq
loci_alt_freq.tbl_df <- function(
    .x,
    # multicore is used by openMP within the
    # freq cpp function
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(nrow(attr(.x$genotypes, "loci")), 1),
    # the bigapply that splits in blocks is not
    # multithreaded, as we use the multiple threads
    # for openMP
    ...) {
  # TODO this is a hack to deal with the class being dropped when going
  # through group_map
  stopifnot_gen_tibble(.x)
  loci_alt_freq(.x$genotypes)
}


#' @export
#' @rdname loci_alt_freq
loci_alt_freq.vctrs_bigSNP <- function(
    .x,
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(nrow(attr(.x, "loci")), 1),
    ...) {
  rlang::check_dots_empty()
  # if we have diploid
  if (is_diploid_only(.x)) {
    loci_alt_freq_diploid(.x, n_cores = n_cores, block_size = block_size)
  } else {
    loci_alt_freq_polyploid(.x, n_cores = n_cores, block_size = block_size, ...)
  }
}

#' @export
#' @rdname loci_alt_freq
loci_alt_freq.grouped_df <- function(
    .x,
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(nrow(attr(.x, "loci")), 1),
    ...) {
  rlang::check_dots_empty()
  if (is_diploid_only(.x)) {
    geno_fbm <- .gt_get_bigsnp(.x)$genotypes
    # rows (individuals) that we want to use
    rows_to_keep <- vctrs::vec_data(.x$genotypes)

    # internal function that can be used with a big_apply
    gt_group_alt_freq_freq_sub <- function(BM, ind, rows_to_keep) {
      freq_mat <- gt_grouped_alt_freq_diploid(
        BM = BM,
        rowInd = rows_to_keep,
        colInd = ind,
        groupIds = dplyr::group_indices(.x) - 1,
        ngroups = max(dplyr::group_indices(.x)),
        ncores = n_cores
      )$freq_alt
    }
    freq_mat <- bigstatsr::big_apply(
      geno_fbm,
      a.FUN = gt_group_alt_freq_freq_sub,
      rows_to_keep = rows_to_keep,
      ind = show_loci(.x)$big_index,
      ncores = 1, # we only use 1 cpu, we let openMP use multiple cores
      # in the cpp code
      block.size = block_size,
      a.combine = "rbind"
    )
    # return a list to mimic a group_map
    lapply(seq_len(ncol(freq_mat)), function(i) freq_mat[, i])
  } else {
    # TODO this is seriously inefficient
    # we should replace it with a cpp function
    group_map(
      .x,
      .f = ~ loci_alt_freq(
        .x, ,
        n_cores = n_cores,
        block_size = block_size,
        ...
      )
    )
  }
}

#' @rdname loci_alt_freq
#' @export
loci_maf <- function(.x, n_cores, block_size, ...) {
  UseMethod("loci_maf", .x)
}

#' @export
#' @rdname loci_alt_freq
loci_maf.tbl_df <- function(
    .x,
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(nrow(attr(.x$genotypes, "loci")), 1),
    ...) {
  # TODO this is a hack to deal with the class being dropped when going
  # through group_map
  stopifnot_gen_tibble(.x)
  loci_maf(.x$genotypes, n_cores = n_cores, block_size = block_size, ...)
}

#' @export
#' @rdname loci_alt_freq
loci_maf.vctrs_bigSNP <- function(
    .x,
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(nrow(attr(.x, "loci")), 1),
    ...) {
  freq <- loci_alt_freq(.x, n_cores = n_cores, block_size = block_size, ...)
  freq[freq > 0.5 & !is.na(freq)] <- 1 - freq[freq > 0.5 & !is.na(freq)]
  freq
}

#' @export
#' @rdname loci_alt_freq
loci_maf.grouped_df <- function(
    .x,
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(nrow(attr(.x, "loci")), 1),
    ...) {
  rlang::check_dots_empty()
  if (is_diploid_only(.x)) {
    geno_fbm <- .gt_get_bigsnp(.x)$genotypes
    # rows (individuals) that we want to use
    rows_to_keep <- vctrs::vec_data(.x$genotypes)

    # internal function that can be used with a big_apply
    gt_group_alt_freq_freq_sub <- function(BM, ind, rows_to_keep) {
      freq_mat <- gt_grouped_alt_freq_diploid(
        BM = BM,
        rowInd = rows_to_keep,
        colInd = ind,
        groupIds = dplyr::group_indices(.x) - 1,
        ngroups = max(dplyr::group_indices(.x)),
        ncores = n_cores
      )$freq_alt
    }
    freq_mat <- bigstatsr::big_apply(
      geno_fbm,
      a.FUN = gt_group_alt_freq_freq_sub,
      rows_to_keep = rows_to_keep,
      ind = show_loci(.x)$big_index,
      ncores = 1, # we only use 1 cpu, we let openMP use multiple cores
      # in the cpp code
      block.size = block_size,
      a.combine = "rbind"
    )
    freq_mat[freq_mat > 0.5 & !is.na(freq_mat)] <-
      1 - freq_mat[freq_mat > 0.5 & !is.na(freq_mat)]
    # return a list to mimic a group_map
    lapply(seq_len(ncol(freq_mat)), function(i) freq_mat[, i])
  } else {
    # the polyploid case
    # TODO this is seriously inefficient
    # we should replace it with a cpp function
    group_map(
      .x,
      .f = ~ loci_maf(.x, n_cores = n_cores, block_size = block_size, ...)
    )
  }
}

# function to estimate frequencies for diploid
loci_alt_freq_diploid <- function(.x, n_cores, block_size) {
  # get the FBM
  geno_fbm <- attr(.x, "bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep) > 1) {
    # create function to use in big_apply
    big_sub_counts <- function(X, ind, rows_to_keep) {
      # nolint
      col_counts <- bigstatsr::big_counts(
        X,
        ind.row = rows_to_keep,
        ind.col = ind
      )
      means_from_counts <- function(x) {
        (x[2] + x[3] * 2) / ((x[1] + x[2] + x[3]) * 2)
      }
      freq_sub <- apply(col_counts, 2, means_from_counts)
      freq_sub
    }
    freq <- bigstatsr::big_apply(
      geno_fbm,
      a.FUN = big_sub_counts,
      rows_to_keep = rows_to_keep,
      ind = attr(.x, "loci")$big_index,
      ncores = 1, # parallelisation is used within the function
      block.size = block_size,
      a.combine = "c"
    )
  } else {
    # if we have a single individual
    freq <- geno_fbm[rows_to_keep, attr(.x, "loci")$big_index] / 2
  }
  freq
}

loci_alt_freq_polyploid <- function(.x, n_cores, block_size, ...) {
  warning(paste(
    "this function still needs a proper unit test!!!",
    "It assumes alleles are the unit of observation"
  ))
  # get the FBM
  geno_fbm <- attr(.x, "bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  ploidy_by_indiv <- indiv_ploidy(.x)
  if (length(rows_to_keep) > 1) {
    # col means for submatrix (all rows, only some columns)
    col_sums_na <- function(X, ind, rows_to_keep, ploidy_by_indiv) {
      # nolint
      res <- colSums(X[rows_to_keep, ind], na.rm = TRUE)
      col_na <- function(a, ploidy_by_indiv) {
        sum(is.na(a) * ploidy_by_indiv)
      }
      res <- cbind(
        res,
        apply(
          X[rows_to_keep, ind],
          2,
          col_na,
          ploidy_by_indiv = ploidy_by_indiv
        )
      )
      res
    }
    col_sums_na_mat <- bigstatsr::big_apply(
      geno_fbm,
      a.FUN = col_sums_na,
      rows_to_keep = rows_to_keep,
      ind = attr(.x, "loci")$big_index,
      ploidy_by_indiv = ploidy_by_indiv,
      ncores = 1, # parallelisation is used within the function
      block.size = block_size,
      a.combine = "rbind"
    )
    # now get frequency accounting for missing values
    col_sums_na_mat[, 1] / (sum(ploidy_by_indiv) - col_sums_na_mat[, 2])
  } else {
    # if we have a single individual
    geno_fbm[rows_to_keep, attr(.x, "loci")$big_index] / ploidy_by_indiv
  }
}
