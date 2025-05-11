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
#' @param type type of object to return One of "tidy" or "pairwise" for a
#'   pairwise matrix of populations. Default is "tidy".
#' @param by_locus_type type of object to return. One of "tidy", "matrix" or
#' "list". Default is "tidy".
#' @param by_locus boolean, determining whether Fst should be returned by
#'   locus(TRUE), or as a single genome wide value obtained by taking the ratio
#'   of the mean numerator and denominator (FALSE, the default).
#' @param method one of 'Hudson', 'Nei87', and 'WC84'
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param return_num_dem returns a list of numerators and denominators for each
#'  locus. This is useful for creating windowed estimates of Fst (as
#'  we need to compute the mean numerator and denominator within
#'  each window). Default is
#'  FALSE. Currently only implemented for Hudson's method.
#' @returns if `type=tidy`, a tibble of genome-wide pairwise Fst values with
#'   each pairwise combination as a row if "by_locus=FALSE", else a list
#'   including the tibble of genome-wide values as well as a matrix with
#'   pairwise Fst by locus with loci as rows and and pairwise combinations as
#'   columns. If `type=pairwise`, a matrix of genome-wide pairwise Fst values is
#'   returned.
#' @export
#' @examples
#' example_gt <- example_gt("gen_tbl")
#'
#' # For a basic global pairwise Fst calculation:
#' example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87")
#'
#' # With a pairwise matrix:
#' example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87", type = "pairwise")
#'
#' # To calculate Fst by locus:
#' example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Hudson", by_locus = TRUE)
#'
pairwise_pop_fst <- function(
    .x,
    type = c("tidy", "pairwise"),
    by_locus = FALSE,
    by_locus_type = c("tidy", "matrix", "list"),
    method = c("Hudson", "Nei87", "WC84"),
    return_num_dem = FALSE,
    n_cores = bigstatsr::nb_cores()) {
  if (n_cores > 1) {
    # Remove checking for two levels of parallelism
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = TRUE))
  }

  type <- match.arg(type)
  if (by_locus) {
    by_locus_type <- match.arg(by_locus_type)
  }

  if (!inherits(.x, "grouped_df")) {
    stop(".x should be a grouped df")
  }
  method <- match.arg(method)

  if (!is.logical(return_num_dem)) {
    stop("return_num_dem must be a logical value (TRUE or FALSE)")
  }

  # if we want to return the numerator and denominator, we need to
  # compute the pairwise Fst for each locus
  if (return_num_dem && !by_locus) {
    message("`by_locus` set to TRUE because `return_num_dem = TRUE`.")
    by_locus <- TRUE
  }

  if (method == "Hudson") {
    pairwise_pop_fst_hudson(
      .x = .x, by_locus = by_locus, type = type, by_locus_type = by_locus_type,
      return_num_dem = return_num_dem, n_cores = n_cores
    )
  } else if (method == "Nei87") {
    pairwise_pop_fst_nei87_cpp(
      .x = .x, by_locus = by_locus,
      type = type, by_locus_type = by_locus_type,
      return_num_dem = return_num_dem, n_cores = n_cores
    )
  } else if (method == "WC84") {
    pairwise_pop_fst_wc84_cpp(
      .x = .x, by_locus = by_locus,
      type = type, by_locus_type = by_locus_type,
      return_num_dem = return_num_dem, n_cores = n_cores
    )
  }
}


pairwise_pop_fst_hudson <- function(
    .x,
    type = type,
    by_locus = FALSE,
    by_locus_type = by_locus_type,
    return_num_dem = FALSE,
    n_cores = bigstatsr::nb_cores()) {
  # get the populations
  .group_levels <- .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels), 2)

  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(
    .gt_get_bigsnp(.x)$genotypes,
    rowInd = .gt_bigsnp_rows(.x),
    colInd = .gt_bigsnp_cols(.x),
    groupIds = dplyr::group_indices(.x) - 1,
    ngroups = nrow(.group_levels),
    ncores = n_cores
  )

  fst_list <- pairwise_fst_hudson_loop(
    pairwise_combn = pairwise_combn,
    n = pop_freqs_df$n,
    freq_alt = pop_freqs_df$freq_alt,
    freq_ref = pop_freqs_df$freq_ref,
    by_locus = by_locus,
    return_num_dem = return_num_dem
  )

  format_fst_list(
    fst_list = fst_list, .x = .x, pairwise_combn = pairwise_combn,
    .group_levels = .group_levels,
    type = type, by_locus_type = by_locus_type, by_locus = by_locus,
    return_num_dem = return_num_dem
  )
}

# the implementation for Nei 87, adapted from hierfstat
pairwise_pop_fst_nei87_cpp <- function(
    .x,
    type = type,
    by_locus = FALSE,
    by_locus_type = by_locus_type,
    return_num_dem = FALSE,
    n_cores = bigstatsr::nb_cores()) {
  # get the populations
  .group_levels <- .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels), 2)
  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(
    .gt_get_bigsnp(.x)$genotypes,
    rowInd = .gt_bigsnp_rows(.x),
    colInd = .gt_bigsnp_cols(.x),
    groupIds = dplyr::group_indices(.x) - 1,
    ngroups = nrow(.group_levels),
    ncores = n_cores
  )
  fst_list <- pairwise_fst_nei87_loop(
    pairwise_combn = pairwise_combn,
    n = pop_freqs_df$n,
    het_obs = pop_freqs_df$het_obs,
    freq_alt = pop_freqs_df$freq_alt,
    freq_ref = pop_freqs_df$freq_ref,
    by_locus = by_locus,
    return_num_dem = return_num_dem
  )
  
  format_fst_list(
    fst_list = fst_list, .x = .x, pairwise_combn = pairwise_combn,
    .group_levels = .group_levels,
    type = type, by_locus_type = by_locus_type, by_locus = by_locus,
    return_num_dem = return_num_dem
  )
}

# the implementation for Nei 87, adapted from hierfstat
pairwise_pop_fst_nei87 <- function(
    .x,
    type = type,
    by_locus = FALSE,
    by_locus_type = by_locus_type,
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

  # format nicely the fst_tot objects
  group_combinations <- cbind(
    .group_levels[pairwise_combn[1, ], ],
    .group_levels[pairwise_combn[2, ], ]
  )
  names(group_combinations) <- c(
    paste0(dplyr::group_vars(.x), "_1"),
    paste0(dplyr::group_vars(.x), "_2")
  )
  fst_tot <- tibble::tibble(group_combinations, value = fst_tot)

  if (type == "pairwise") { # if we return a matrix
    fst_tot <- tidy_to_matrix(fst_tot)
  }
  if (by_locus && by_locus_type == "matrix") {
    rownames(fst_locus) <- loci_names(.x)
    colnames(fst_locus) <- col_names_combinations(group_combinations,
      prefix = "fst"
    )
    return(list(Fst_by_locus = fst_locus, Fst = fst_tot))
  } else if (by_locus && by_locus_type == "tidy") {
    fst_mat_tbl <- as.data.frame(fst_locus)
    colnames(fst_mat_tbl) <- col_names_combinations(group_combinations,
      prefix = "fst"
    )
    fst_mat_tbl$loci <- loci_names(.x)
    cols <- names(fst_mat_tbl)[names(fst_mat_tbl) != "loci"]
    long_fst_loc <- fst_mat_tbl %>%
      tidyr::pivot_longer(cols = all_of(cols), names_to = "stat_name")
    return(list(Fst_by_locus = long_fst_loc, Fst = fst_tot))
  } else if (by_locus && by_locus_type == "list") {
    fst_mat_tbl <- as.data.frame(fst_locus)
    colnames(fst_mat_tbl) <- col_names_combinations(group_combinations,
      prefix = "fst"
    )
    fst_locus <- as.list(fst_mat_tbl)
    return(list(Fst_by_locus = fst_locus, Fst = fst_tot))
  } else {
    return(fst_tot)
  }
}


# This should be equivalent to the hierfstat implementation
pairwise_pop_fst_wc84 <- function(
    .x,
    type = type,
    by_locus = FALSE,
    by_locus_type = by_locus_type,
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

  if (type == "pairwise") { # if we return a matrix
    fst_tot <- tidy_to_matrix(fst_tot)
  }
  if (by_locus && by_locus_type == "matrix") {
    rownames(fst_locus) <- loci_names(.x)
    colnames(fst_locus) <- col_names_combinations(group_combinations,
      prefix = "fst"
    )
    return(list(Fst_by_locus = fst_locus, Fst = fst_tot))
  } else if (by_locus && by_locus_type == "tidy") {
    fst_mat_tbl <- as.data.frame(fst_locus)
    colnames(fst_mat_tbl) <- col_names_combinations(group_combinations,
      prefix = "fst"
    )
    fst_mat_tbl$loci <- loci_names(.x)
    cols <- names(fst_mat_tbl)[names(fst_mat_tbl) != "loci"]
    long_fst_loc <- fst_mat_tbl %>%
      tidyr::pivot_longer(cols = all_of(cols), names_to = "stat_name")
    return(list(Fst_by_locus = long_fst_loc, Fst = fst_tot))
  } else if (by_locus && by_locus_type == "list") {
    fst_mat_tbl <- as.data.frame(fst_locus)
    colnames(fst_mat_tbl) <- col_names_combinations(group_combinations,
      prefix = "fst"
    )
    fst_locus <- as.list(fst_mat_tbl)
    return(list(Fst_by_locus = fst_locus, Fst = fst_tot))
  } else {
    return(fst_tot)
  }
}

# This should be equivalent to the hierfstat implementation
pairwise_pop_fst_wc84_cpp <- function(
    .x,
    type = type,
    by_locus = FALSE,
    by_locus_type,
    return_num_dem = FALSE,
    n_cores = bigstatsr::nb_cores()) {
  # get the populations
  .group_levels <- .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels), 2)

  # summarise population frequencies
  pop_freqs_df <- gt_grouped_summaries(
    .gt_get_bigsnp(.x)$genotypes,
    rowInd = .gt_bigsnp_rows(.x),
    colInd = .gt_bigsnp_cols(.x),
    groupIds = dplyr::group_indices(.x) - 1,
    ngroups = nrow(.group_levels),
    ncores = n_cores
  )
  fst_list <- pairwise_fst_wc84_loop(
    pairwise_combn = pairwise_combn,
    n = pop_freqs_df$n,
    freq_alt = pop_freqs_df$freq_alt,
    het_obs = pop_freqs_df$het_obs,
    by_locus = by_locus,
    return_num_dem = return_num_dem
  )
  format_fst_list(
    fst_list, .x, pairwise_combn, .group_levels,
    type, by_locus_type, by_locus
  )
}


# This function is used to turn a tidy tibble into a matrix
tidy_to_matrix <- function(tidy_tbl) {
  fst_tot_wide <- tidyr::pivot_wider(
    tidy_tbl,
    names_from = "population_2",
    values_from = "value"
  ) %>%
    tibble::column_to_rownames(var = "population_1") %>%
    as.matrix()
  # add missing row and col to make the matrix symmetrical
  fst_tot_wide <- cbind(NA_real_, fst_tot_wide)
  fst_tot_wide <- rbind(fst_tot_wide, NA_real_)
  # fix dim names
  rownames(fst_tot_wide)[nrow(fst_tot_wide)] <-
    utils::tail(colnames(fst_tot_wide), 1)
  colnames(fst_tot_wide)[1] <- rownames(fst_tot_wide)[1]
  # fill lower triangle
  fst_tot_wide[lower.tri(fst_tot_wide)] <-
    t(fst_tot_wide)[lower.tri(fst_tot_wide)]
  return(fst_tot_wide)
}

# This function is used to create the column names for the pairwise
# combinations
col_names_combinations <- function(group_combinations, prefix = NULL) {
  comb_labels <- apply(
    group_combinations,
    1,
    function(x) paste(x, collapse = ".")
  )
  if (!is.null(prefix)) {
    comb_labels <- paste(prefix, comb_labels, sep = "_")
  }
  return(comb_labels)
}

# Format output for pairwise_fst_ functions based on type

format_fst_list <- function(fst_list, .x, pairwise_combn, .group_levels,
                            type, by_locus_type, by_locus, return_num_dem = FALSE) {
  # format nicely the group combinations object
  # (we'll use it to create column names)
  group_combinations <- cbind(
    .group_levels[pairwise_combn[1, ], ],
    .group_levels[pairwise_combn[2, ], ]
  )
  names(group_combinations) <- c(
    paste0(dplyr::group_vars(.x), "_1"),
    paste0(dplyr::group_vars(.x), "_2")
  )
  # if we want the numerator and denominator,
  # we need to format them and return them
  if (return_num_dem) {
    rownames(fst_list$Fst_by_locus_num) <- loci_names(.x)
    colnames(fst_list$Fst_by_locus_num) <-
      col_names_combinations(group_combinations)
    rownames(fst_list$Fst_by_locus_den) <- loci_names(.x)
    colnames(fst_list$Fst_by_locus_den) <-
      col_names_combinations(group_combinations)
    return(fst_list)
  }

  # otherwise we start formatting the other objects
  fst_tot <- tibble::tibble(group_combinations, value = fst_list$fst_tot)

  if (type == "pairwise") { # if we return a matrix
    fst_tot <- tidy_to_matrix(fst_tot)
  }

  if (by_locus && by_locus_type == "matrix") {
    rownames(fst_list$fst_locus) <- loci_names(.x)
    colnames(fst_list$fst_locus) <- col_names_combinations(group_combinations,
      prefix = "fst"
    )
    return(list(Fst_by_locus = fst_list$fst_locus, Fst = fst_tot))
  } else if (by_locus && by_locus_type == "tidy") {
    fst_mat_tbl <- as.data.frame(fst_list$fst_locus)
    colnames(fst_mat_tbl) <- col_names_combinations(group_combinations,
      prefix = "fst"
    )
    fst_mat_tbl$loci <- loci_names(.x)
    cols <- names(fst_mat_tbl)[names(fst_mat_tbl) != "loci"]
    long_fst_loc <- fst_mat_tbl %>%
      tidyr::pivot_longer(cols = all_of(cols), names_to = "stat_name")
    return(list(Fst_by_locus = long_fst_loc, Fst = fst_tot))
  } else if (by_locus && by_locus_type == "list") {
    fst_mat_tbl <- as.data.frame(fst_list$fst_locus)
    colnames(fst_mat_tbl) <- col_names_combinations(group_combinations,
      prefix = "fst"
    )
    fst_list$fst_locus <- as.list(fst_mat_tbl)
    return(list(Fst_by_locus = fst_list$fst_locus, Fst = fst_tot))
  } else {
    return(fst_tot)
  }
}
