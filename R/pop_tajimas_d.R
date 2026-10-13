#' Estimate Tajima's D for the whole genome
#'
#' Note that Tajima's D estimates from data that have been filtered or
#' ascertained can be difficult to interpret. This function should ideally
#' be used on sequence data prior to filtering.
#'
#' Nucleotide diversity is estimated at each locus with [loci_pi()], so this
#' function works with diploid, polyploid, mixed-ploidy and pseudohaploid data.
#' The sample size is the total number of allele copies sampled (i.e. the sum
#' of the ploidy of the individuals, with pseudohaploids contributing one
#' copy each). Loci with missing nucleotide diversity (e.g. loci with no
#' genotyped individuals) are ignored.
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
    # get pi for each locus (this deals with polyploids and pseudohaploids)
    pi <- loci_pi(.x, n_cores = n_cores, block_size = block_size)
    # number of sampled allele copies (assuming all individuals contribute
    # some data)
    n <- sum(indiv_ploidy(.x))
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
  ...
) {
  rlang::check_dots_empty()
  # get pi for each locus and group (this deals with polyploids and
  # pseudohaploids)
  pi_by_grp <- loci_pi(.x,
    n_cores = n_cores, block_size = block_size,
    type = "list"
  )
  # number of sampled allele copies per group
  n_by_grp <- .n_allele_copies(.x)
  # compute tajimas d for each group
  tajimas_d <- Map(tajimas_d_from_pi_vec, unname(pi_by_grp), n_by_grp)
  # return a list to mimic a group_map
  return(tajimas_d)
}


# estimate tajimas d from a vector of pi per locus and the number of
# sampled allele copies (i.e. 2 * N individuals for diploids). Loci with
# missing pi are ignored
tajimas_d_from_pi_vec <- function(pi, n) {
  # count segregating sites from pi
  seg <- sum(pi > 0.0 & pi < 1.0, na.rm = TRUE)
  k_hat <- sum(pi, na.rm = TRUE)
  return(tajimas_d_from_sums(k_hat = k_hat, seg = seg, n = n))
}

# estimate tajimas d from the sum of pi across loci (k_hat), the number of
# segregating sites (seg) and the number of haploid genomes (n). This is
# vectorised over k_hat and seg (e.g. to compute leave-one-out estimates
# for a block jackknife)
tajimas_d_from_sums <- function(k_hat, seg, n) {
  a1 <- sum(1 / (1:(n - 1)))
  a2 <- sum(1 / ((1:(n - 1))^2))
  e1 <- ((n + 1) / (3 * (n - 1)) - 1 / a1) / a1
  e2_num <- 2 * (n^2 + n + 3) / (9 * n * (n - 1)) -
    (n + 2) / (n * a1) + a2 / (a1^2)
  e2_den <- a1^2 + a2
  vd <- e1 * seg + e2_num / e2_den * seg * (seg - 1)
  return((k_hat - seg / a1) / sqrt(vd))
}

# Get the number of sampled allele copies for each group (or for the whole
# gen_tibble if it is not grouped). Pseudohaploid genotypes are stored as 0/2
# dosages but contribute one allele per sample.
.n_allele_copies <- function(.x) {
  if (inherits(.x, "grouped_gen_tbl")) {
    .x %>%
      dplyr::summarise(n = sum(indiv_ploidy(.data$genotypes))) %>%
      dplyr::pull(.data$n)
  } else {
    sum(indiv_ploidy(.x))
  }
}
