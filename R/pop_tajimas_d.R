#' Estimate Tajima's D for the whole genome
#'
#' Note that Tajima's D estimates from data that have been filtered or
#' ascertained can be difficult to interpret. This function should ideally
#' be used on sequence data prior to filtering.
#'
#' Nucleotide diversity is estimated at each locus with [loci_pi()], so this
#' function works with diploid, polyploid, mixed-ploidy and pseudohaploid data.
#'
#' Missing genotypes are accounted for following Bailey et al. (2025), as
#' implemented in `pixy`: the sample size \eqn{n_x} at each locus \eqn{x} is
#' the number of non-missing allele copies at that locus (i.e. the sum of the
#' ploidy of the genotyped individuals, with pseudohaploids contributing one
#' copy each). The two estimators of \eqn{\theta} are then:
#' \itemize{
#'   \item \eqn{\hat{\theta}_\pi = \sum_x \pi_x}, where \eqn{\pi_x} is
#'     the nucleotide diversity at locus \eqn{x}, estimated from its
#'     \eqn{n_x} allele copies (loci with fewer than 2 allele copies
#'     contribute 0);
#'   \item \eqn{\hat{\theta}_W = \sum_x 1 / a_{n_x}}, summed over the
#'     segregating sites, where \eqn{a_{n_x} = \sum_{i=1}^{n_x - 1} 1/i} (i.e.
#'     each segregating site is weighted by its own sample size).
#' }
#' The standard deviation of \eqn{\hat{\theta}_\pi - \hat{\theta}_W} is
#' computed with the formula of Tajima (1989), using the number of
#' segregating sites \eqn{S} and the mean sample size across loci with at
#' least one non-missing allele copy (rounded to the nearest integer, as the
#' formula is only defined for integer sample sizes). When there are no
#' missing data, all loci have the same sample size and the statistic is
#' identical to the classic Tajima's D. Loci with no genotyped individuals
#' are ignored. Tajima's D is `NA` when there are no segregating sites, or
#' when the mean sample size is smaller than 2.
#'
#' @references Tajima F (1989) Statistical method for testing the neutral
#'   mutation hypothesis by DNA polymorphism. Genetics 123:585-595.
#'
#'   Bailey N, Stevison L, Samuk K (2025) Correcting for bias in estimates of
#'   \eqn{\theta_W} and Tajima's D from missing data in next-generation
#'   sequencing. Molecular Ecology Resources, e14104.
#'   \doi{10.1111/1755-0998.14104}
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
  ...
) {
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
  ...
) {
  rlang::check_dots_empty()

  # as long as we have more than one individual
  if (length(.x) > 1) {
    # get pi and the number of valid allele copies for each locus (this deals
    # with polyploids and pseudohaploids), treating all individuals as a
    # single group
    pi_n <- .loci_pi_n(.x,
      group_ids = rep(0, length(.x)),
      n_groups = 1,
      n_cores = n_cores,
      block_size = block_size
    )
    tajimas_d <- tajimas_d_from_pi_vec(pi_n$pi[, 1], pi_n$n[, 1])
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
  ...
) {
  rlang::check_dots_empty()
  # get pi and the number of valid allele copies for each locus and group
  # (this deals with polyploids and pseudohaploids)
  pi_n <- .loci_pi_n(.x$genotypes,
    group_ids = dplyr::group_indices(.x) - 1,
    n_groups = max(dplyr::group_indices(.x)),
    n_cores = n_cores,
    block_size = block_size
  )
  # compute tajimas d for each group
  tajimas_d <- lapply(seq_len(ncol(pi_n$pi)), function(i_grp) {
    tajimas_d_from_pi_vec(pi_n$pi[, i_grp], pi_n$n[, i_grp])
  })
  # return a list to mimic a group_map
  return(tajimas_d)
}


# estimate tajimas d from a vector of pi per locus and the number of
# valid allele copies used to estimate pi at each locus (n). n can be a single
# value (e.g. 2 * N individuals for diploids with no missing data), or a
# vector with one value per locus. Loci with no valid allele copies are ignored
tajimas_d_from_pi_vec <- function(pi, n) {
  # get the quantities needed to compute tajimas d for each locus
  locus_terms <- tajimas_d_locus_terms(pi, n)
  # sum them across loci, and compute tajimas d
  tajimas_d_from_term_sums(colSums(locus_terms, na.rm = TRUE))
}

# estimate tajimas d from the sums of the quantities computed by
# tajimas_d_locus_terms(); term_sums is either a named vector (one estimate),
# or a matrix with named columns (one estimate per row, e.g. for each window
# or jackknife block)
tajimas_d_from_term_sums <- function(term_sums) {
  # recast a vector as a matrix with a single row
  if (is.null(dim(term_sums))) {
    term_sums <- matrix(term_sums,
      nrow = 1,
      dimnames = list(NULL, names(term_sums))
    )
  }
  tajimas_d <- tajimas_d_from_sums(
    k_hat = term_sums[, "k_hat"],
    seg = term_sums[, "seg"],
    theta_w = term_sums[, "theta_w"],
    n_sum = term_sums[, "n"],
    n_loci = term_sums[, "n_loci"]
  )
  return(unname(tajimas_d))
}

# compute, for each locus, the quantities that are summed across loci to
# estimate tajimas d (following Bailey et al. 2025, each locus has its own
# sample size n, i.e. the number of valid allele copies):
# - k_hat: pi at the locus (0 if there are fewer than 2 valid allele copies)
# - seg: 1 if the locus is segregating, 0 otherwise
# - theta_w: the contribution of the locus to Watterson's theta, i.e.
#   1 / sum_{i=1}^{n-1} 1/i for a segregating site, 0 otherwise
# - n: the sample size of the locus
# - n_loci: 1 for a locus with at least one valid allele copy
# All values are NA for loci with no valid allele copies, so that they can be
# ignored when summing (e.g. across a window or a jackknife block)
tajimas_d_locus_terms <- function(pi, n) {
  # recycle n if a single value was given
  n <- rep_len(n, length(pi))
  # loci with no valid allele copies are not used
  observed <- !is.na(n) & n > 0
  # pi is not defined with fewer than 2 allele copies; such loci are
  # not segregating, and contribute 0 to the sum of pi (as in pixy)
  pi[observed & is.na(pi)] <- 0
  # pi is greater than 0 only if both alleles are present
  seg <- as.numeric(pi > 0.0)
  # contribution to watterson's theta (only segregating sites contribute)
  theta_w <- ifelse(seg > 0, 1 / harmonic_sum(n, power = 1), 0)
  cbind(
    k_hat = ifelse(observed, pi, NA),
    seg = ifelse(observed, seg, NA),
    theta_w = ifelse(observed, theta_w, NA),
    n = ifelse(observed, n, NA),
    n_loci = ifelse(observed, 1, NA)
  )
}

# estimate tajimas d from the sum across loci of the quantities computed by
# tajimas_d_locus_terms():
# - k_hat: the sum of pi
# - seg: the number of segregating sites
# - theta_w: watterson's theta, i.e. the sum over segregating sites of
#   1 / sum_{i=1}^{n-1} 1/i (Bailey et al. 2025, equation 3)
# - n_sum: the sum of the sample size of each locus
# - n_loci: the number of loci
# This is vectorised over all arguments (e.g. to compute estimates for
# windows, or leave-one-out estimates for a block jackknife). With no missing
# data (i.e. the same n for all loci), this is the classic Tajima's D
tajimas_d_from_sums <- function(k_hat, seg, theta_w, n_sum, n_loci) {
  # the standard deviation of Tajima (1989) is computed for the mean sample
  # size across loci (rounded to the nearest integer, as the formula is only
  # defined for integer sample sizes), as done in pixy
  n <- round(n_sum / n_loci)
  # the formula is not defined for fewer than 2 allele copies
  n[n < 2] <- NA
  a1 <- harmonic_sum(n, power = 1)
  a2 <- harmonic_sum(n, power = 2)
  e1 <- ((n + 1) / (3 * (n - 1)) - 1 / a1) / a1
  e2_num <- 2 * (n^2 + n + 3) / (9 * n * (n - 1)) -
    (n + 2) / (n * a1) + a2 / (a1^2)
  e2_den <- a1^2 + a2
  vd <- e1 * seg + e2_num / e2_den * seg * (seg - 1)
  tajimas_d <- (k_hat - theta_w) / sqrt(vd)
  # tajimas d is not defined without segregating sites
  tajimas_d[seg == 0] <- NA
  return(tajimas_d)
}

# compute sum_{i=1}^{n-1} 1/i^power for each value of n (a vector of
# integers, possibly with NAs); this is a1 (power = 1) and a2 (power = 2) in
# Tajima (1989)
harmonic_sum <- function(n, power = 1) {
  res <- rep(NA_real_, length(n))
  # only compute the sums for valid values of n (the sum is 0 for n = 1)
  valid <- !is.na(n) & n >= 1
  if (any(valid)) {
    # cumulative sums up to the largest n
    n_max <- max(n[valid])
    cum_sums <- c(0, cumsum(1 / (seq_len(n_max - 1))^power))
    # index the cumulative sums by n
    res[valid] <- cum_sums[n[valid]]
  }
  res
}

# get pi and the number of valid allele copies for each locus and each group
# of a gen_tibble (if the gen_tibble is not grouped, all individuals are
# treated as a single group). Returns a list with two matrices (pi and n), with
# one row per locus and one column per group
.pi_n_by_group <- function(.x,
                           n_cores = bigstatsr::nb_cores(),
                           block_size = bigstatsr::block_size(nrow(.x), 1)) {
  if (inherits(.x, "grouped_df")) {
    group_ids <- dplyr::group_indices(.x) - 1
  } else {
    group_ids <- rep(0, nrow(.x))
  }
  .loci_pi_n(.x$genotypes,
    group_ids = group_ids,
    n_groups = max(group_ids) + 1,
    n_cores = n_cores,
    block_size = block_size
  )
}
