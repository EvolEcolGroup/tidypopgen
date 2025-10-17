#' Estimate Tajima's D for the whole genome
#'
#' Note that Tajima's D estimates from data that have been filtered or
#' ascertained can be difficult to interpret. This function should ideally
#' be used on sequence data prior to filtering.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param block_size maximum number of loci read at once.
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns A single numeric value (Tajima's D D) for the whole data set, `NA`
#'   when the statistic is not defined. For grouped data a list of Tajima's D D
#'   values (one per group) is returned.
#' @rdname pop_tajimas_d
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # Compute Tajima's D
#' example_gt %>% pop_tajimas_d()
pop_tajimas_d <- function(.x, n_cores, block_size, ...) {
  UseMethod("pop_tajimas_d", .x)
}

#' @export
#' @rdname pop_tajimas_d
pop_tajimas_d.tbl_df <- function(
    .x,
    # multicore is used by openMP within the
    # freq cpp function
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(nrow(.x), 1),
    # the bigapply that splits in blocks is not
    # multithreaded, as we use the multiple threads
    # for openMP
    ...) {
  # TODO this is a hack to deal with the class being dropped when going
  # through group_map
  stopifnot_gen_tibble(.x)
  pop_tajimas_d(.x$genotypes, n_cores = n_cores, block_size = block_size)
}


#' @export
#' @rdname pop_tajimas_d
pop_tajimas_d.vctrs_bigSNP <- function(
    .x,
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(length(.x), 1),
    ...) {
  rlang::check_dots_empty()

  stopifnot_diploid(.x)
  # if we have diploid
  # get the FBM
  geno_fbm <- .gt_get_fbm(.x)
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep) > 1) {
    # internal function that can be used with a big_apply #nolint start
    gt_pi_sub <- function(BM, ind, rows_to_keep) {
      gt_pi_diploid(
        BM = BM,
        rowInd = rows_to_keep,
        colInd = ind,
        ncores = n_cores
      )
    } # nolint end
    pi <- bigstatsr::big_apply(
      geno_fbm,
      a.FUN = gt_pi_sub,
      rows_to_keep = rows_to_keep,
      ind = attr(.x, "loci")$big_index,
      ncores = 1, # parallelisation is used within the function
      block.size = block_size,
      a.combine = "c"
    )
    # n individuals (assuming they all contribute some data)
    n <- length(.x) * 2
    tajimas_d <- tajimas_d_from_pi_vec(pi, n)
  } else {
    # if we have a single individual
    # tajimas_d does not really make sense for a single individual
    tajimas_d <- NA
  }
  tajimas_d
}

#' @export
#' @rdname pop_tajimas_d
pop_tajimas_d.grouped_df <- function(
    .x,
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(nrow(.x), 1),
    ...) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  geno_fbm <- .gt_get_fbm(.x)
  # rows (individuals) that we want to use
  rows_to_keep <- .gt_fbm_rows(.x)

  # internal function that can be used with a big_apply #nolint start
  gt_group_pi_sub <- function(BM, ind, rows_to_keep) {
    freq_mat <- gt_grouped_pi_diploid(
      BM = BM,
      rowInd = rows_to_keep,
      colInd = ind,
      groupIds = dplyr::group_indices(.x) - 1,
      ngroups = max(dplyr::group_indices(.x)),
      ncores = n_cores
    )$pi
  } # nolint end
  pi_mat <- bigstatsr::big_apply(
    geno_fbm,
    a.FUN = gt_group_pi_sub,
    rows_to_keep = rows_to_keep,
    ind = show_loci(.x)$big_index,
    ncores = 1, # we only use 1 cpu, we let openMP use multiple cores
    # in the cpp code
    block.size = block_size,
    a.combine = "rbind"
  )
  # number of rows per group
  n_by_grp <- .x %>%
    dplyr::summarise(n = n()) %>%
    dplyr::pull(.data$n)
  tajimas_d <- list()
  # TODO this should be done with an apply function or in cpp
  for (i_grp in seq_along(n_by_grp)) {
    # get the number of individuals in the group
    n <- n_by_grp[i_grp] * 2
    # compute tajimas d for each group
    tajimas_d[i_grp] <- tajimas_d_from_pi_vec(pi_mat[, i_grp], n)
  }
  # return a list to mimic a group_map
  return(tajimas_d)
}


# estimate tajimas d from a vector of pi per locus and the number of
# haploid genomes (i.e. 2 * N individuals)
tajimas_d_from_pi_vec <- function(pi, n) {
  # count segregating sites from pi
  seg <- sum(pi > 0.0 & pi < 1.0, na.rm = TRUE)
  k_hat <- sum(pi)

  a1 <- sum(1 / (1:(n - 1)))
  a2 <- sum(1 / ((1:(n - 1))^2))
  e1 <- ((n + 1) / (3 * (n - 1)) - 1 / a1) / a1
  e2_num <- 2 * (n^2 + n + 3) / (9 * n * (n - 1)) -
    (n + 2) / (n * a1) + a2 / (a1^2)
  e2_den <- a1^2 + a2
  vd <- e1 * seg + e2_num / e2_den * seg * (seg - 1)
  return((k_hat - seg / a1) / sqrt(vd))
}
