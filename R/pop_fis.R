#' Compute population specific FIS
#'
#' This function computes population specific FIS, using either the approach of
#' Nei 1987 (with an algorithm equivalent to the one used by
#' `hierfstat::basic.stats()`) or of Weir and Goudet 2017 (with an algorithm
#' equivalent to the one used by `hierfstat::fis.dosage()`).
#' @references Nei M. (1987) Molecular Evolutionary Genetics. Columbia
#'   University Press.
#' @references Weir, BS and Goudet J (2017) A Unified Characterization of
#'   Population Structure and Relatedness. Genetics (2017) 206:2085
#' @param .x a grouped [`gen_tibble`] (as obtained by using
#'   [dplyr::group_by()]); if method == "Nei87", it is also possible to use
#'   a regular `gen_tibble`, in which case all individuals are considered to
#'   belong to a single population.
#' @param method one of "Nei87" (based on Nei 1987, eqn 7.41) or "WG17" (for
#'   Weir and Goudet 2017) to compute FIS
#' @param by_locus boolean, determining whether FIS should be returned by
#'   locus(TRUE), or as a single genome wide value (FALSE, the default). Note
#'   that this is only relevant for "Nei87", as "WG17" always returns a single
#'   value.
#' @param include_global boolean determining whether, besides the population
#'   specific estimates, a global estimate should be appended. Note that this
#'   will return a vector of n populations plus 1 (the global value), or a
#'   matrix with n+1 columns if `by_locus=TRUE`.
#' @param allele_sharing_mat optional and only relevant for "WG17", the matrix
#'   of Allele Sharing returned by [pairwise_allele_sharing()] with
#'   `as_matrix=TRUE`. As a number of statistics can be derived from the Allele
#'   Sharing matrix, it is sometimes more efficient to pre-compute this matrix.
#' @param n_cores number of cores to be used for "Nei87" (it defaults to
#'   [bigstatsr::nb_cores()]); "WG17" uses the allele sharing matrix, see
#'   [pairwise_allele_sharing()].
#' @returns a vector of population specific fis (plus the global value if
#'   `include_global=TRUE`)
#' @export
#' @seealso [hierfstat::basic.stats()] [hierfstat::fis.dosage()]
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # Compute FIS using Nei87:
#' example_gt %>% pop_fis(method = "Nei87")
#'
#' # Compute FIS using WG17:
#' example_gt %>% pop_fis(method = "WG17")
#'
#' # To include the global FIS, set include_global = TRUE:
#' example_gt %>% pop_fis(method = "Nei87", include_global = TRUE)
#'
#' # To return FIS by locus, set by_locus = TRUE:
#' example_gt %>% pop_fis(method = "Nei87", by_locus = TRUE)
#'
#' # To calculate from a pre-computed allele sharing matrix:
#' allele_sharing_mat <- pairwise_allele_sharing(example_gt, as_matrix = TRUE)
#' example_gt %>% pop_fis(
#'   method = "WG17", allele_sharing_mat =
#'     allele_sharing_mat
#' )
pop_fis <- function(
  .x,
  method = c("Nei87", "WG17"),
  by_locus = FALSE,
  include_global = FALSE,
  allele_sharing_mat = NULL,
  n_cores = bigstatsr::nb_cores()
) {
  method <- match.arg(method)
  if (method == "Nei87") {
    if (!is.null(allele_sharing_mat)) {
      stop("allele_sharing_mat not relevant for Nei87")
    }
    pop_fis_nei87(.x,
      by_locus = by_locus, include_global = include_global,
      n_cores = n_cores
    )
  } else if (method == "WG17") {
    if (by_locus) {
      stop("by_locus not implemented for WG17")
    }
    pop_fis_wg17(
      .x,
      include_global = include_global,
      allele_sharing_mat = allele_sharing_mat
    )
  }
}

pop_fis_nei87 <- function(
  .x,
  by_locus = FALSE,
  include_global = FALSE,
  n_cores = bigstatsr::nb_cores()
) {
  stopifnot_diploid(.x)
  # compute Fis = 1 - Ho/Hs for each population (and the global statistics,
  # which include the global Fis) in a single pass over the genotypes
  pop_stats <- .pop_stats_dip_pseudo(
    .x,
    by_locus_fis = by_locus,
    by_locus_global = by_locus && include_global,
    n_cores = n_cores
  )

  if (by_locus) {
    Fis <- pop_stats$fis # nolint
    colnames(Fis) <- pop_stats$group_names # nolint
    if (include_global) {
      # Fis is the ninth global statistic
      global <- pop_stats$global[, 9]
      Fis <- cbind(Fis, global) # nolint
    }
  } else {
    # means across loci
    Fis <- pop_stats$mean_fis # nolint
    names(Fis) <- pop_stats$group_names
    if (include_global) {
      # as in pop_global_stats(), the genome-wide global Fis is computed from
      # the means of Ho and Hs (ignoring infinite values)
      global <- 1 - pop_stats$mean_global_finite[1] /
        pop_stats$mean_global_finite[2]
      names(global) <- "global"
      Fis <- c(Fis, global) # nolint
    }
  }
  return(Fis)
}


pop_fis_wg17 <- function(
  .x,
  include_global = FALSE,
  allele_sharing_mat = NULL
) {
  stopifnot_diploid(.x)
  if (!inherits(.x, "grouped_df")) {
    stop(".x should be a grouped gen_tibble")
  }
  if (is.null(allele_sharing_mat)) {
    allele_sharing_mat <- pairwise_allele_sharing(.x, as_matrix = TRUE)
  }

  Mij <- allele_sharing_mat
  Mii <- diag(Mij) * 2 - 1
  diag(Mij) <- NA
  pop <- factor(group_indices(.x))
  pop_levels <- levels(pop)
  n_pop <- length(pop_levels)
  if (n_pop < 2) {
    stop("pop_fis with method 'WG17' requires at least two groups")
  }
  wil <- lapply(pop_levels, function(z) which(pop == z))
  Fi <- lapply(wil, function(pop_levels) Mii[pop_levels])
  Fsts <- unlist(lapply(
    wil,
    function(pop_levels) {
      mean(Mij[pop_levels, pop_levels], na.rm = TRUE)
    }
  ))
  Mb <- 0
  mMij <- matrix(numeric(n_pop^2), ncol = n_pop)
  for (i in 2:n_pop) {
    p1 <- wil[[i]]
    for (j in 1:(i - 1)) {
      p2 <- wil[[j]]
      mMij[i, j] <- mMij[j, i] <- mean(Mij[p1, p2], na.rm = TRUE)
      Mb <- Mb + mMij[i, j]
    }
  }
  diag(mMij) <- Fsts

  Mb <- Mb * 2 / (n_pop * (n_pop - 1))

  # estimate individual Fi (inbreeding)
  for (i in 1:n_pop) {
    Fi[[i]] <- (Fi[[i]] - Fsts[[i]]) / (1 - Fsts[[i]])
  }
  names(Fi) <- as.character(pop_levels)
  fis_by_pop <- unlist(lapply(Fi, mean, na.rm = TRUE))

  if (include_global) {
    fis_by_pop <- c(fis_by_pop, mean(fis_by_pop, na.rm = TRUE))
    names(fis_by_pop) <- c(dplyr::group_keys(.x) %>% pull(1), "global")
  } else {
    names(fis_by_pop) <- c(dplyr::group_keys(.x) %>% pull(1))
  }
  return(fis_by_pop)
}
