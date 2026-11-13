#' Compute the Identity by State Matrix for a bigSNP object
#'
#' This function computes the IBS matrix.
#'
#' Note that monomorphic sites are counted (as in PLINK). Remove monomorphic
#' sites beforehand if this is a concern.
#'
#' See [pairwise_ibs()] for the two corrections for missing genotypes
#' (`missing_correction`).
#' @param X a [bigstatsr::FBM.code256] matrix (as found in the `genotypes`
#' slot of a [bigsnpr::bigSNP] object). Genotypes are decoded with its
#' `code256`, and values other than 0, 1 and 2 are treated as missing.
#' @param ind.row An optional vector of the row indices that are used. If not
#'   specified, all rows are used. Don't use negative indices.
#' @param ind.col An optional vector of the column indices that are used. If not
#'   specified, all columns are used. Don't use negative indices.
#' @param type one of:
#'   - "proportion" (the default): the proportion of alleles that are identical
#'   by state, a similarity (as `--distance ibs flat-missing` in PLINK, which
#'   writes it to a .mibs file).
#'   - "adjusted_counts": the allele count distance, a dissimilarity (as
#'   `--distance flat-missing` in PLINK, which writes it to a .dist file): the
#'   number of alleles that differ between two individuals, scaled up by the
#'   number of loci over the number of loci genotyped in both individuals, to
#'   correct for missing genotypes. It equals `2 * n_loci * (1 -
#'   proportion)`, where `n_loci` is the number of loci used, so it is 0 for
#'   identical individuals and at most `2 * n_loci`. Note that, without
#'   `flat-missing`, PLINK corrects for missing genotypes differently, which
#'   gives slightly different values when genotypes are missing (see
#'   `missing_correction`).
#'   - "raw_counts": the counts of identical alleles and of non-missing
#'   alleles, from which the two other quantities are computed.
#' @param block.size maximum number of columns read at once. Genotypes are
#'   compared in a compact bit-packed form, so this has little effect on
#'   memory use.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param missing_correction the correction for missing genotypes of
#'   "proportion" and "adjusted_counts": "flat" (the default, as PLINK 1.9's
#'   `flat-missing` modifier) or "allele_freq" (PLINK 1.9's default, which
#'   weights the missing loci by their allele frequencies). See
#'   [pairwise_ibs()] for details. It must be "flat" for "raw_counts".
#' @param allele_freqs only used if `missing_correction = "allele_freq"`: the
#'   frequency of the alternate allele at each locus in `ind.col` (as PLINK's
#'   `--read-freq`); either allele can be used, as the weights are the same
#'   for both. If `NULL` (the default), it is computed from the genotypes of
#'   the individuals in `ind.row`. Missing values are set to 0.5, as PLINK
#'   does for loci without genotypes.
#' @returns for "proportion" and "adjusted_counts", a matrix with one row and
#'   one column per individual; for "raw_counts", a list of two
#'   [bigstatsr::FBM] matrices, `ibs` (the counts of alleles identical by
#'   state) and `valid_n` (the counts of non-missing alleles, i.e. twice the
#'   number of loci genotyped in both individuals).
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' X <- attr(example_gt$genotypes, "fbm")
#' snp_ibs(X)
#'
#' # Compute for individuals 1 to 5
#' snp_ibs(X, ind.row = 1:5, ind.col = 1:5)
#'
#' # Adjust block.size
#' snp_ibs(X, block.size = 2)
#'
#' # Change type
#' snp_ibs(X, type = "proportion")
#' snp_ibs(X, type = "adjusted_counts")
#' snp_ibs(X, type = "raw_counts")
#'
#' # PLINK 1.9's default correction for missing genotypes
#' snp_ibs(X, type = "adjusted_counts", missing_correction = "allele_freq")
#'
snp_ibs <- function(
  X, # nolint start
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  type = c("proportion", "adjusted_counts", "raw_counts"),
  block.size = bigstatsr::block_size(nrow(X)),
  n_cores = bigstatsr::nb_cores(),
  missing_correction = c("flat", "allele_freq"),
  allele_freqs = NULL
) {
  # nolint end
  type <- match.arg(type)
  missing_correction <- match.arg(missing_correction)
  weighted <- missing_correction == "allele_freq"
  if (weighted && type == "raw_counts") {
    stop(
      "missing_correction = 'allele_freq' can not be used with ",
      "type = 'raw_counts'"
    )
  }
  if (!is.null(allele_freqs) && !weighted) {
    stop("allele_freqs is only used with missing_correction = 'allele_freq'")
  }

  n <- length(ind.row)
  # no individuals (e.g. an empty group in group_map()): an empty matrix, as
  # FBMs can not be empty
  if (n == 0) {
    if (type == "raw_counts") {
      stop("type = 'raw_counts' needs at least one individual")
    }
    return(matrix(numeric(0), 0, 0))
  }
  # the count FBMs are temporary unless they are returned (raw_counts): make
  # them in a private directory, deleted on exit, so that their names cannot
  # clash with other .bk files (see with_private_fbm_dir())
  if (type != "raw_counts") {
    local_private_fbm_dir()
  }
  # FBM matrix to count the IBS counts
  IBS <- bigstatsr::FBM(n, n, init = 0) # nolint
  # FBM matrix to store valid number of comparisons (i.e NOT NA)
  IBS_valid_loci <- bigstatsr::FBM(n, n, init = 0) # nolint
  m <- length(ind.col)

  if (weighted) {
    # PLINK 1.9's integer weights of missing loci (see plink_missing_weights())
    if (is.null(allele_freqs)) {
      allele_freqs <- snp_alt_freq(X, ind.row = ind.row, ind.col = ind.col)
    }
    missing_wts <- plink_missing_weights(allele_freqs, m)
    # FBM of the sums of the weights of the loci missing in both individuals
    IBS_missing_wts <- bigstatsr::FBM(n, n, init = 0) # nolint
  }

  intervals <- CutBySize(m, block.size)

  # for each block of loci, genotypes are packed into bits and compared
  # with popcounts in C++ (see src/pairwise_bitpacked.h), so no dense
  # copy of the block is needed
  for (j in bigstatsr::rows_along(intervals)) {
    ind <- seq2(intervals[j, ]) # this iteration indices
    ind.col.ind <- ind.col[ind] # subset given indices by the iteration indices #nolint
    if (weighted) {
      increment_ibs_counts_weighted(
        IBS,
        IBS_valid_loci,
        IBS_missing_wts,
        X,
        ind.row,
        ind.col.ind,
        wts = missing_wts[ind],
        ncores = n_cores
      )
    } else {
      increment_ibs_counts(
        IBS,
        IBS_valid_loci,
        X,
        ind.row,
        ind.col.ind,
        ncores = n_cores
      )
    }
  }

  if (weighted) {
    # PLINK 1.9's weighted correction: the number of alleles that differ,
    # scaled up by the total weight of all loci over the total weight of the
    # loci genotyped in both individuals. The operations are those of PLINK
    # (calc_distance() in plink_calc.c), so the values are identical to
    # PLINK's if the weights are.
    w_total <- sum(missing_wts)
    w_missing <- IBS_missing_wts[]
    w_missing_ind <- diag(w_missing)
    # the weight of the loci genotyped in both individuals:
    # all - missing in i - missing in j + missing in both
    w_valid <- w_total - outer(w_missing_ind, w_missing_ind, "+") + w_missing
    # number of alleles that differ (valid alleles - identical alleles)
    n_differ <- IBS_valid_loci[] - IBS[]
    ibs_dist <- (w_total / w_valid) * n_differ
    if (type == "adjusted_counts") {
      return(ibs_dist)
    } else {
      # as PLINK's --distance ibs: 1 - the distance as a proportion of the
      # 2 * m alleles
      return(1 - ibs_dist * (0.5 / m))
    }
  }

  if (type == "raw_counts") {
    return(list(ibs = IBS, valid_n = IBS_valid_loci))
  } else {
    # IBS proportion
    # get the means of each column
    divide_sub <- function(X, ind, Y) (X[, ind] / Y[, ind]) # nolint
    ibs_prop <- bigstatsr::big_apply(
      IBS,
      a.FUN = divide_sub,
      Y = IBS_valid_loci,
      a.combine = "cbind"
    )
    if (type == "proportion") {
      return(ibs_prop)
    } else {
      # the allele count distance (as PLINK --distance flat-missing): the
      # number of alleles that differ, scaled up to all m loci to correct for
      # missing genotypes, i.e. 2 * m * (1 - ibs_prop)
      ibs_dist <- 2 * m * (1 - ibs_prop)
      return(ibs_dist)
    }
  }
}

#' PLINK 1.9's weights of missing loci
#'
#' The weights of the loci used by PLINK 1.9 to correct distances for
#' missing genotypes (without `flat-missing`, see [pairwise_ibs()]): the
#' expected distance between two individuals at the locus under Hardy-Weinberg
#' equilibrium, up to a constant, i.e. `p * (1 - p) * (p^2 - p + 1)` for an
#' allele frequency `p` (0 for monomorphic loci), scaled to sum to just under
#' `2^32` (`2^32 - n_loci`) and rounded to integers, as in `calc_distance()`
#' in `plink_calc.c` of PLINK 1.9.
#' @param allele_freqs the frequency of either allele at each locus (missing
#'   values are set to 0.5, as PLINK does for loci without genotypes)
#' @param n_loci the number of loci
#' @returns a vector of integer weights (stored as doubles, as they can be
#'   larger than the largest R integer)
#' @keywords internal
#' @noRd
plink_missing_weights <- function(allele_freqs, n_loci) {
  if (length(allele_freqs) != n_loci) {
    stop("allele_freqs should have one value per locus")
  }
  if (any(allele_freqs < 0 | allele_freqs > 1, na.rm = TRUE)) {
    stop("allele_freqs should be between 0 and 1")
  }
  p <- allele_freqs
  p[is.na(p)] <- 0.5
  # expected distance under HWE, up to a constant (exactly 0 if monomorphic)
  raw_wts <- ifelse(p == 0 | p == 1, 0, p * (1.0 - p) * (p * p - p + 1))
  raw_sum <- sum(raw_wts)
  if (!(raw_sum > 0)) {
    stop(
      "all loci are monomorphic, so the weights of missing_correction = ",
      "'allele_freq' are all zero"
    )
  }
  # scale to sum to just under 2^32 (minus n_loci, to guard against rounding
  # up) and round to integers, as PLINK's (uint32_t)(x + 0.5)
  floor(raw_wts * ((4294967296.0 - n_loci) / raw_sum) + 0.5)
}

#' Frequency of the alternate allele from an FBM.code256
#'
#' The frequency of the alternate allele (the allele counted by the
#' genotypes) at each locus, from the genotypes decoded with the code256 of
#' the FBM (so imputed genotypes are used only if the code256 decodes them);
#' values other than 0, 1 and 2 are treated as missing, as in snp_ibs().
#' @param X an FBM.code256
#' @param ind.row,ind.col the rows (individuals) and columns (loci) used
#' @returns a vector of frequencies (NaN for loci without genotypes)
#' @keywords internal
#' @noRd
snp_alt_freq <- function(X, ind.row, ind.col) { # nolint
  counts <- bigstatsr::big_counts(X, ind.row = ind.row, ind.col = ind.col)
  # rows of counts are the decoded values (the last one, NA, has no name)
  values <- suppressWarnings(as.numeric(rownames(counts)))
  count_of <- function(g) {
    if (g %in% values) counts[which(values == g), ] else rep(0, ncol(counts))
  }
  n_1 <- count_of(1)
  n_2 <- count_of(2)
  n_valid <- count_of(0) + n_1 + n_2
  (n_1 + 2 * n_2) / (2 * n_valid)
}
