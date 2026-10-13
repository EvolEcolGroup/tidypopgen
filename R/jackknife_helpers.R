# Internal helpers for the weighted block jackknife (Busing et al. 1999)

#' Assign loci to jackknife blocks
#'
#' Blocks are non-overlapping and never span chromosomes. Block ids are
#' integers that are unique across the whole genome.
#'
#' @param loci_table a dataframe including at least a column 'chromosome', and
#'   additionally a column 'position' if `size_unit` is "bp", or a column
#'   'genetic_dist' if `size_unit` is "cM".
#' @param block_size the size of the blocks
#' @param size_unit one of "snp", "bp" or "cM"
#' @returns an integer vector of block ids, one per locus
#' @keywords internal
#' @noRd
.jackknife_block_ids <- function(loci_table, block_size,
                                 size_unit = c("snp", "bp", "cM")) {
  size_unit <- match.arg(size_unit)
  .check_loci_table_units(loci_table, size_unit)
  # check that block_size is positive
  if (!is.numeric(block_size) || length(block_size) != 1 || block_size <= 0) {
    stop("block_size must be a single positive number.")
  }
  # integer index for each chromosome (we avoid creating character vectors
  # with one element per locus, as they are very large for big datasets)
  chrom_idx <- match(loci_table$chromosome, unique(loci_table$chromosome))
  # block id within each chromosome
  if (size_unit == "snp") {
    # count the snps within each chromosome
    snp_idx <- stats::ave(seq_along(chrom_idx), chrom_idx, FUN = seq_along)
    local_id <- ceiling(snp_idx / block_size)
  } else if (size_unit == "bp") {
    local_id <- floor(loci_table$position / block_size)
  } else if (size_unit == "cM") {
    local_id <- floor(loci_table$genetic_dist / block_size)
  }
  # combine chromosome and local block id into a single numeric key (local
  # ids are shifted to start from 0, so that keys are unique across
  # chromosomes), and number the blocks in order of appearance
  local_id <- local_id - min(local_id)
  block_key <- chrom_idx * (max(local_id) + 1) + local_id
  match(block_key, unique(block_key))
}

#' Sum a matrix of per locus values within each block
#'
#' @param mat a matrix (or vector) of per locus values (loci as rows)
#' @param block_id a vector of block ids (as returned by
#'   `.jackknife_block_ids()`)
#' @returns a list with `sums` (a matrix with one row per block, with the sum
#'   of the non-NA values in each block) and `n` (a matrix with the number of
#'   non-NA values in each block)
#' @keywords internal
#' @noRd
.block_sums <- function(mat, block_id) {
  # remove all attributes other than dim (e.g. names of the loci, or
  # matrix_int_names), only keeping the column names (one per statistic); we
  # avoid as.matrix(), which would copy the names of a vector as row names
  mat <- matrix(as.numeric(mat),
    nrow = NROW(mat), ncol = NCOL(mat),
    dimnames = list(NULL, colnames(mat))
  )
  sums <- rowsum(mat, block_id, na.rm = TRUE, reorder = TRUE)
  n <- rowsum((!is.na(mat)) * 1, block_id, reorder = TRUE)
  # rowsum() names the rows with the block ids (as characters), which we drop
  # (rows are in the order of the block ids)
  rownames(sums) <- NULL
  rownames(n) <- NULL
  list(sums = sums, n = n)
}

#' Weighted block jackknife
#'
#' Compute the weighted delete-m jackknife estimate and standard error, as
#' described in Busing et al. (1999). This function is vectorised over
#' statistics (columns).
#'
#' @param theta_full a vector of estimates based on all the data (one per
#'   statistic)
#' @param theta_loo a matrix of leave-one-out estimates (one row per block,
#'   one column per statistic)
#' @param m a matrix of the number of loci in each block (same dimensions as
#'   `theta_loo`). Blocks with 0 loci are ignored.
#' @returns a tibble with one row per statistic, with columns `estimate`,
#'   `jack_estimate`, `se`, `z` and `n_blocks`
#' @keywords internal
#' @noRd
.jackknife_weighted <- function(theta_full, theta_loo, m) {
  theta_loo <- as.matrix(theta_loo)
  m <- as.matrix(m)
  if (!all(dim(theta_loo) == dim(m))) {
    stop("theta_loo and m must have the same dimensions")
  }
  if (length(theta_full) != ncol(theta_loo)) {
    stop("theta_full must have one value per column of theta_loo")
  }
  jack_estimate <- se <- rep(NA_real_, length(theta_full))
  n_blocks <- integer(length(theta_full))
  for (i_stat in seq_along(theta_full)) {
    # only use blocks that contain data for this statistic
    in_use <- m[, i_stat] > 0
    m_j <- m[in_use, i_stat]
    theta_j <- theta_loo[in_use, i_stat]
    g <- length(m_j)
    n_blocks[i_stat] <- g
    # we need at least 2 blocks, and well defined estimates
    if (g < 2 || !is.finite(theta_full[i_stat]) || !all(is.finite(theta_j))) {
      next
    }
    n <- sum(m_j)
    h_j <- n / m_j
    # jackknife estimate (equation 8 in Busing et al. 1999)
    jack_estimate[i_stat] <- g * theta_full[i_stat] -
      sum((1 - m_j / n) * theta_j)
    # pseudovalues and variance (equations 9 and 10)
    tau_j <- h_j * theta_full[i_stat] - (h_j - 1) * theta_j
    se[i_stat] <- sqrt(
      sum((tau_j - jack_estimate[i_stat])^2 / (h_j - 1)) / g
    )
  }
  tibble::tibble(
    estimate = as.numeric(theta_full),
    jack_estimate = jack_estimate,
    se = se,
    z = jack_estimate / se,
    n_blocks = n_blocks
  )
}

#' Compute full and leave-one-out pairwise Fst for jackknife blocks
#'
#' @param .x a grouped `gen_tibble`
#' @param method the method to use for calculating Fst
#' @param block_id a vector of block ids, one per locus
#' @returns a list with `full` (a named vector of Fst based on all loci),
#'   `loo` (a matrix of leave-one-out Fst, one row per block and one column per
#'   pair of populations, named `fst_a.b`), and `m` (a matrix of the number of
#'   loci used in each block)
#' @keywords internal
#' @noRd
.jackknife_fst_blocks <- function(.x, method, block_id) {
  # create the pairwise Fst by locus, saving numerator and denominator
  pair_fst <- pairwise_pop_fst(.x,
    method = method, return_num_dem = TRUE,
    by_locus = TRUE
  )
  stat_names <- paste0("fst_", colnames(pair_fst$Fst_by_locus_num))
  # convert to plain numeric matrices (dropping the matrix_int_names class)
  num <- matrix(as.numeric(pair_fst$Fst_by_locus_num),
    nrow = nrow(pair_fst$Fst_by_locus_num)
  )
  den <- matrix(as.numeric(pair_fst$Fst_by_locus_den),
    nrow = nrow(pair_fst$Fst_by_locus_den)
  )
  # only use loci where both numerator and denominator are available
  # (as done when computing the genome-wide Fst in pairwise_pop_fst)
  valid <- !is.na(num) & !is.na(den)
  num[!valid] <- NA
  den[!valid] <- NA
  num_blocks <- .block_sums(num, block_id)
  den_blocks <- .block_sums(den, block_id)
  # Fst on all loci, and leaving each block out in turn
  num_tot <- colSums(num_blocks$sums)
  den_tot <- colSums(den_blocks$sums)
  full <- num_tot / den_tot
  loo <- (rep_rows(num_tot, nrow(num_blocks$sums)) - num_blocks$sums) /
    (rep_rows(den_tot, nrow(den_blocks$sums)) - den_blocks$sums)
  names(full) <- stat_names
  colnames(loo) <- stat_names
  m <- num_blocks$n
  colnames(m) <- stat_names
  list(full = full, loo = loo, m = m)
}

# create a matrix with n_rows copies of the vector x as its rows
rep_rows <- function(x, n_rows) {
  matrix(x, nrow = n_rows, ncol = length(x), byrow = TRUE)
}

#' Format jackknife results to mirror the windows functions
#'
#' @param jack a tibble as returned by `.jackknife_weighted()`
#' @param stat_names the names of the statistics
#' @param type either "matrix" (one row per quantity, one column per
#'   statistic) or "tidy" (one row per statistic)
#' @param name_col the name of the column holding the statistic names in the
#'   tidy format
#' @keywords internal
#' @noRd
.jackknife_format <- function(jack, stat_names, type = c("matrix", "tidy"),
                              name_col = "stat_name") {
  type <- match.arg(type)
  if (type == "tidy") {
    res <- dplyr::bind_cols(
      tibble::tibble(!!name_col := stat_names),
      jack
    )
  } else {
    quantities <- names(jack)
    res <- as.data.frame(t(as.matrix(jack)))
    names(res) <- stat_names
    res <- dplyr::bind_cols(
      tibble::tibble(quantity = quantities),
      tibble::as_tibble(res)
    )
  }
  res
}
