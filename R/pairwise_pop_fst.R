#' Compute pairwise population Fst
#'
#' This function computes pairwise Fst. The following methods are implemented:
#' - 'Hudson': Hudson's formulation, as derived in Bhatia et al (2013)
#' for diploids.
#' - 'Nei87' : Fst according to Nei (1987) - includes the correction for
#' heterozygosity when computing Ht, and is equivalent to
#' `hierfstat::pairwise.neifst()`,
#' - 'WC84' : Weir and Cockerham (1984), correcting for missing data and is
#' equivalent to `hierfstat::pairwise.WCfst()`.
#'
#' For all formulae, the genome wide estimate is obtained by taking the ratio of
#' the mean numerators and denominators over all relevant SNPs.
#'
#' @references Bhatia G, Patterson N, Sankararaman S, Price AL. (2013)
#'   Estimating and Interpreting FST: The Impact of Rare Variants. Genome
#'   Research, 23(9):1514–1521.
#'
#'   Nei, M. (1987) Molecular Evolutionary Genetics. Columbia University Press
#'
#'   Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the
#'   analysis of population structure. Evolution, 38(6): 1358–1370.
#'
#' @param .x a grouped [`gen_tibble`] (as obtained by using [dplyr::group_by()])
#' @param tidy boolean whether to return a tidy tibble. Default is TRUE, FALSE
#'   returns a matrix. Set to TRUE if by_locus is TRUE.
#' @param by_locus boolean, determining whether Fst should be returned by
#'   locus(TRUE), or as a single genome wide value obtained by taking the ratio
#'   of the mean numerator and denominator (FALSE, the default).
#' @param method one of 'Hudson', 'Nei87', and 'WC84'
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @returns a tibble of genome-wide pairwise Fst values with each pairwise
#'   combination as a row if "by_locus=FALSE", else a list including the tibble
#'   of genome-wide values as well as a matrix with pairwise Fst by locus with
#'   loci as rows and and pairwise combinations as columns.
#' @export

pairwise_pop_fst <- function(
    .x,
    tidy = TRUE,
    by_locus = FALSE,
    method = c("Hudson", "Nei87", "WC84"),
    n_cores = bigstatsr::nb_cores()) {
  if (n_cores > 1) {
    # Remove checking for two levels of parallelism
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = TRUE))
  }

  if (!inherits(.x, "grouped_df")) {
    stop(".x should be a grouped df")
  }
  if (by_locus) {
    tidy <- TRUE
  }
  method <- match.arg(method)
  if (method == "Hudson") {
    pairwise_pop_fst_hudson(.x = .x, by_locus = by_locus, tidy = tidy)
  } else if (method == "Nei87") {
    pairwise_pop_fst_nei87(.x = .x, by_locus = by_locus, tidy = tidy)
  } else if (method == "WC84") {
    pairwise_pop_fst_wc84(.x = .x, by_locus = by_locus, tidy = tidy)
  }
}

pairwise_pop_fst_hudson <- function(
    .x,
    tidy = TRUE,
    by_locus = FALSE,
    n_cores = bigstatsr::nb_cores()) {
  # get the populations
  .group_levels <- .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels), 2)
  # vector and matrix to store Fst for total and by locus
  fst_tot <- rep(NA_real_, ncol(pairwise_combn))
  if (by_locus) {
    fst_locus <-
      matrix(NA_real_, nrow = count_loci(.x), ncol = ncol(pairwise_combn))
  }
  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(
    .gt_get_bigsnp(.x)$genotypes,
    rowInd = .gt_bigsnp_rows(.x),
    colInd = .gt_bigsnp_cols(.x),
    groupIds = dplyr::group_indices(.x) - 1,
    ngroups = nrow(.group_levels),
    ncores = n_cores
  )

  for (i_col in seq_len(ncol(pairwise_combn))) {
    pop1 <- pairwise_combn[1, i_col]
    pop2 <- pairwise_combn[2, i_col]

    numerator <- (pop_freqs_df$freq_alt[, pop1] -
      pop_freqs_df$freq_alt[, pop2])^2 - # nolint start
      (pop_freqs_df$freq_alt[, pop1] * pop_freqs_df$freq_ref[, pop1]) /
        (pop_freqs_df$n[, pop1] - 1) -
      (pop_freqs_df$freq_alt[, pop2] * pop_freqs_df$freq_ref[, pop2]) /
        (pop_freqs_df$n[, pop2] - 1)
    denominator <- pop_freqs_df$freq_alt[, pop1] *
      pop_freqs_df$freq_ref[, pop2] +
      pop_freqs_df$freq_alt[, pop2] * pop_freqs_df$freq_ref[, pop1] # nolint end
    if (by_locus) {
      fst_locus[, i_col] <- numerator / denominator
    }
    fst_tot[i_col] <- mean(numerator, na.rm = TRUE) /
      mean(denominator, na.rm = TRUE) # nolint
  }
  # format nicely the objects
  group_combinations <- cbind(
    .group_levels[pairwise_combn[1, ], ],
    .group_levels[pairwise_combn[2, ], ]
  )
  names(group_combinations) <- c(
    paste0(dplyr::group_vars(.x), "_1"),
    paste0(dplyr::group_vars(.x), "_2")
  )
  fst_tot <- tibble::tibble(group_combinations, value = fst_tot)

  if (!tidy) {
    fst_tot_wide <- tidyr::spread(fst_tot,
      key = .data$population_2,
      value = .data$value
    )
    matrix_rownames <- c(
      as.vector(fst_tot_wide[1, 1])[[1]],
      names(fst_tot_wide)[-1]
    )
    fst_tot_wide <- fst_tot_wide[, -1]
    fst_tot_wide <- cbind(NA_real_, fst_tot_wide)
    fst_tot_wide <- rbind(fst_tot_wide, NA_real_)
    fst_tot_matrix <- as.matrix(fst_tot_wide)
    rownames(fst_tot_matrix) <- colnames(fst_tot_matrix) <- matrix_rownames
    # if values are present only in lower half matrix, transpose them
    missing_upper <- is.na(fst_tot_matrix) & upper.tri(fst_tot_matrix)
    fst_tot_matrix[missing_upper] <- t(fst_tot_matrix)[missing_upper]
    # then transpose to make the matrix symmetrical
    fst_tot_matrix[lower.tri(fst_tot_matrix)] <-
      t(fst_tot_matrix)[lower.tri(fst_tot_matrix)]
    return(fst_tot_matrix)
  } else if (tidy) {
    if (by_locus) {
      rownames(fst_locus) <- loci_names(.x)
      colnames(fst_locus) <- apply(
        group_combinations,
        1,
        function(x) paste(x, collapse = ".")
      )
      return(list(Fst_by_locus = fst_locus, Fst = fst_tot))
    } else {
      return(fst_tot)
    }
  }
}

## use tidyr::pivot_wider to turn into a matrix if that's what is requested.

# the implementation for Nei 87, adapted from hierfstat
pairwise_pop_fst_nei87 <- function(
    .x,
    tidy = TRUE,
    by_locus = FALSE,
    n_cores = bigstatsr::nb_cores()) {
  # get the populations
  .group_levels <- .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels), 2)
  # vector and matrix to store Fst for total and by locus
  fst_tot <- rep(NA_real_, ncol(pairwise_combn))
  if (by_locus) {
    fst_locus <- matrix(
      NA_real_,
      nrow = count_loci(.x),
      ncol = ncol(pairwise_combn)
    )
  }
  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(
    .gt_get_bigsnp(.x)$genotypes,
    rowInd = .gt_bigsnp_rows(.x),
    colInd = .gt_bigsnp_cols(.x),
    groupIds = dplyr::group_indices(.x) - 1,
    ngroups = nrow(.group_levels),
    ncores = n_cores
  )
  for (i_col in seq_len(ncol(pairwise_combn))) {
    pop1 <- pairwise_combn[1, i_col]
    pop2 <- pairwise_combn[2, i_col]

    n <- pop_freqs_df$n[, c(pop1, pop2)] / 2
    sHo <- pop_freqs_df$het_obs[, c(pop1, pop2)] # nolint
    mHo <- apply(sHo, 1, mean, na.rm = TRUE) # nolint
    freq_alt <- pop_freqs_df$freq_alt[, c(pop1, pop2)]
    freq_ref <- pop_freqs_df$freq_ref[, c(pop1, pop2)]

    # sum of squared frequencies
    sp2 <- freq_alt^2 + freq_ref^2
    Hs <- (1 - sp2 - sHo / 2 / n) # nolint
    Hs <- n / (n - 1) * Hs # nolint
    np <- apply(n, 1, fun <- function(x) sum(!is.na(x))) # nolint
    # mean sample size over the populations
    mn <- apply(
      n,
      1,
      fun <- function(x) {
        sum(!is.na(x)) / sum(1 / x[!is.na(x)])
      }
    )
    # mean sum of square frequencies
    msp2 <- apply(sp2, 1, mean, na.rm = TRUE) # nolint start
    mp2 <- rowMeans(freq_alt)^2 + rowMeans(freq_ref)^2
    mHs <- mn / (mn - 1) * (1 - msp2 - mHo / 2 / mn)
    Ht <- 1 - mp2 + mHs / mn / np - mHo / 2 / mn / np

    Dst <- Ht - mHs
    Dstp <- np / (np - 1) * Dst
    Htp <- mHs + Dstp
    if (by_locus) {
      fst_locus[, i_col] <- Dstp / Htp
    }
    fst_tot[i_col] <- mean(Dstp) / mean(Htp) # nolint end
  }
  # format nicely the objects
  group_combinations <- cbind(
    .group_levels[pairwise_combn[1, ], ],
    .group_levels[pairwise_combn[2, ], ]
  )
  names(group_combinations) <- c(
    paste0(dplyr::group_vars(.x), "_1"),
    paste0(dplyr::group_vars(.x), "_2")
  )
  fst_tot <- tibble::tibble(group_combinations, value = fst_tot)
  if (!tidy) {
    fst_tot_wide <- tidyr::spread(fst_tot,
      key = .data$population_2,
      value = .data$value
    )
    matrix_rownames <- c(
      as.vector(fst_tot_wide[1, 1])[[1]],
      names(fst_tot_wide)[-1]
    )
    fst_tot_wide <- fst_tot_wide[, -1]
    fst_tot_wide <- cbind(NA_real_, fst_tot_wide)
    fst_tot_wide <- rbind(fst_tot_wide, NA_real_)
    fst_tot_matrix <- as.matrix(fst_tot_wide)
    rownames(fst_tot_matrix) <- colnames(fst_tot_matrix) <- matrix_rownames
    # if values are present only in lower half matrix, transpose them
    missing_upper <- is.na(fst_tot_matrix) & upper.tri(fst_tot_matrix)
    fst_tot_matrix[missing_upper] <- t(fst_tot_matrix)[missing_upper]
    # then transpose to make the matrix symmetrical
    fst_tot_matrix[lower.tri(fst_tot_matrix)] <-
      t(fst_tot_matrix)[lower.tri(fst_tot_matrix)]
    return(fst_tot_matrix)
  } else if (tidy) {
    if (by_locus) {
      rownames(fst_locus) <- loci_names(.x)
      colnames(fst_locus) <- apply(
        group_combinations,
        1,
        function(x) paste(x, collapse = ".")
      )
      return(list(Fst_by_locus = fst_locus, Fst = fst_tot))
    } else {
      return(fst_tot)
    }
  }
}


# This should be equivalent to the hierfstat implementation
pairwise_pop_fst_wc84 <- function(
    .x,
    tidy = TRUE,
    by_locus = FALSE,
    n_cores = bigstatsr::nb_cores()) {
  # get the populations
  .group_levels <- .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels), 2)
  # vector and matrix to store Fst for total and by locus
  fst_tot <- rep(NA_real_, ncol(pairwise_combn))
  if (by_locus) {
    fst_locus <- matrix(
      NA_real_,
      nrow = count_loci(.x),
      ncol = ncol(pairwise_combn)
    )
  }
  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(
    .gt_get_bigsnp(.x)$genotypes,
    rowInd = .gt_bigsnp_rows(.x),
    colInd = .gt_bigsnp_cols(.x),
    groupIds = dplyr::group_indices(.x) - 1,
    ngroups = nrow(.group_levels),
    ncores = n_cores
  )
  for (i_col in seq_len(ncol(pairwise_combn))) {
    pops <- pairwise_combn[c(1, 2), i_col]
    # allele counts for this pair of populations
    r <- length(pops) # number of populations (to generalise the function)
    an <- pop_freqs_df$n[, pops]
    # number of individuals
    n_ind <- an / 2
    n_total <- rowSums(n_ind)
    n_bar <- rowMeans(n_ind)
    n_c <- (n_total - (rowSums(n_ind^2) / n_total)) / (r - 1)

    # allele frequencies for this pair of populations
    p <- pop_freqs_df$freq_alt[, pops]
    # mean frequencies
    p_bar <- rowSums(p * n_ind) / n_total
    s2 <- rowSums((p - p_bar)^2 * n_ind) / n_bar / (r - 1)
    h_bar <- rowSums(pop_freqs_df$het_obs[, pops] * n_ind) / n_total
    # nolint start
    a <- n_bar /
      n_c *
      (s2 -
        1 / (n_bar - 1) * (p_bar * (1 - p_bar) - (r - 1) / r * s2 - h_bar / 4))
    b <- n_bar /
      (n_bar - 1) *
      (p_bar *
        (1 - p_bar) -
        (r - 1) / r * s2 -
        (2 * n_bar - 1) / (4 * n_bar) * h_bar) # nolint
    c <- h_bar / 2
    # nolint end
    numerator <- a
    denominator <- a + b + c
    if (by_locus) {
      fst_locus[, i_col] <- numerator / denominator
    }
    fst_tot[i_col] <-
      mean(numerator, na.rm = TRUE) / mean(denominator, na.rm = TRUE)
  }
  # format nicely the objects
  group_combinations <- cbind(
    .group_levels[pairwise_combn[1, ], ],
    .group_levels[pairwise_combn[2, ], ]
  )
  names(group_combinations) <- c(
    paste0(dplyr::group_vars(.x), "_1"),
    paste0(dplyr::group_vars(.x), "_2")
  )
  fst_tot <- tibble::tibble(group_combinations, value = fst_tot)

  if (!tidy) {
    fst_tot_wide <- tidyr::spread(fst_tot,
      key = .data$population_2,
      value = .data$value
    )
    matrix_rownames <- c(
      as.vector(fst_tot_wide[1, 1])[[1]],
      names(fst_tot_wide)[-1]
    )
    fst_tot_wide <- fst_tot_wide[, -1]
    fst_tot_wide <- cbind(NA_real_, fst_tot_wide)
    fst_tot_wide <- rbind(fst_tot_wide, NA_real_)
    fst_tot_matrix <- as.matrix(fst_tot_wide)
    rownames(fst_tot_matrix) <- colnames(fst_tot_matrix) <- matrix_rownames
    # if values are present only in lower half matrix, transpose them
    missing_upper <- is.na(fst_tot_matrix) & upper.tri(fst_tot_matrix)
    fst_tot_matrix[missing_upper] <- t(fst_tot_matrix)[missing_upper]
    # then transpose to make the matrix symmetrical
    fst_tot_matrix[lower.tri(fst_tot_matrix)] <-
      t(fst_tot_matrix)[lower.tri(fst_tot_matrix)]
    return(fst_tot_matrix)
  } else if (tidy) {
    if (by_locus) {
      rownames(fst_locus) <- loci_names(.x)
      colnames(fst_locus) <- apply(
        group_combinations,
        1,
        function(x) paste(x, collapse = ".")
      )
      return(list(Fst_by_locus = fst_locus, Fst = fst_tot))
    } else {
      return(fst_tot)
    }
  }
}
