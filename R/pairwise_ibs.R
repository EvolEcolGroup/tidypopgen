#' Compute the Identity by State Matrix for a `gen_tibble` object
#'
#' This function computes the IBS matrix.
#'
#' Note that monomorphic sites are currently considered (as in PLINK). Remove
#' monomorphic sites before running `pairwise_ibs()` if this is a concern.
#'
#' Imputed genotypes are used if `x` is set to use them (see
#' [gt_set_imputed()]), and treated as missing otherwise.
#'
#' Each pair of individuals is compared only on the loci genotyped in both.
#' Two corrections for missing genotypes are available (`missing_correction`),
#' which give the same values when no genotypes are missing:
#' - "flat" (the default), as PLINK 1.9 with the `flat-missing` modifier
#' (`--distance ibs flat-missing` for "proportion", `--distance
#' flat-missing` for "adjusted_counts"): all loci have the same weight, so
#' the number of alleles that differ between individuals `i` and `j`,
#' `d_ij`, is scaled up by `n_loci / n_ij`, where `n_ij` is the number of loci
#' genotyped in both ("adjusted_counts"), and "proportion" is
#' `1 - d_ij / (2 * n_ij)`.
#' - "allele_freq", as PLINK 1.9's default (`--distance ibs` and
#' `--distance`): loci are weighted by their expected contribution to the
#' distance, which depends on their allele frequency `p`. Under Hardy-Weinberg
#' equilibrium, the expected number of alleles that differ between two
#' individuals at a locus is `4 * p * (1 - p) * (p^2 - p + 1)`, and PLINK
#' uses `w = p * (1 - p) * (p^2 - p + 1)` as the weight of the locus (0 for
#' monomorphic loci). "adjusted_counts" is then `d_ij * W / W_ij`, where `W`
#' is the sum of the weights of all loci and `W_ij` that of the loci genotyped
#' in both individuals, and "proportion" is `1 - adjusted_counts / (2 *
#' n_loci)`. As in PLINK, the weights are scaled to sum to just under `2^32`
#' and rounded to integers (so values are identical to PLINK's), `W_ij`
#' excludes each locus missing in either individual or in both (its weight is
#' removed once), and monomorphic loci have no weight (a missing genotype at a
#' monomorphic locus does not change the distance, but the locus still counts
#' in `n_loci`). PLINK 1.9 computes the allele frequencies from the founders
#' only (unless `--nonfounders` is used), or reads them with `--read-freq`;
#' `pairwise_ibs()` computes them from all individuals of `x` (the same as
#' PLINK when all individuals are founders, as in data without a pedigree),
#' unless they are given with `allele_freqs` (e.g. from the founders, with
#' [loci_alt_freq()]). Loci without any genotype have a frequency of 0.5, as
#' in PLINK. PLINK's `1-ibs` modifier (`--distance 1-ibs`, written to a `.mdist`
#' file) gives `1 - proportion` with either correction.
#'
#' Pairs of individuals without any locus genotyped in both (or, with
#' "allele_freq", without any polymorphic locus genotyped in both) have
#' `NaN` (or `Inf`) values.
#' @param x a `gen_tibble` object.
#' @param as_matrix boolean, determining whether the results should be a square
#'   symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)
#' @param type one of:
#'   - "proportion" (the default): the proportion of alleles that are identical
#'   by state, a similarity (as `--distance ibs flat-missing` in PLINK, which
#'   writes it to a .mibs file; with `missing_correction = "allele_freq"`, as
#'   `--distance ibs`, see Details).
#'   - "adjusted_counts": the allele count distance, a dissimilarity (as
#'   `--distance flat-missing` in PLINK, which writes it to a .dist file): the
#'   number of alleles that differ between two individuals, scaled up by the
#'   number of loci over the number of loci genotyped in both individuals, to
#'   correct for missing genotypes. It equals `2 * n_loci * (1 -
#'   proportion)`, so it is 0 for identical individuals and at most `2 *
#'   n_loci`. With `missing_correction = "allele_freq"`, as `--distance` in
#'   PLINK (without `flat-missing`), see Details.
#'   - "raw_counts": the counts of identical alleles and of non-missing
#'   alleles, from which the two other quantities are computed.
#' @param block_size maximum number of loci read at once. Genotypes are
#'   compared in a compact bit-packed form, so this has little effect on
#'   memory use.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param missing_correction the correction for missing genotypes, for
#'   "proportion" and "adjusted_counts": "flat" (the default, as PLINK 1.9
#'   with the `flat-missing` modifier) or "allele_freq" (as PLINK 1.9's
#'   default, which weights loci by their allele frequencies); see Details.
#'   It must be "flat" for "raw_counts".
#' @param allele_freqs only used if `missing_correction = "allele_freq"`: the
#'   frequency of the alternate allele at each locus (as PLINK's
#'   `--read-freq`). If `NULL` (the default), it is computed from all the
#'   individuals of `x` (using imputed genotypes if `x` is set to use them).
#' @returns for "proportion" and "adjusted_counts", a tidy tibble with one row
#'   per pair of individuals (columns `item1`, `item2` and `value`), or a
#'   square matrix if `as_matrix = TRUE`; for "raw_counts", a list of two
#'   [bigstatsr::FBM] matrices, `ibs` (the counts of alleles identical by
#'   state) and `valid_n` (the counts of non-missing alleles, i.e. twice the
#'   number of loci genotyped in both individuals).
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' pairwise_ibs(example_gt, type = "proportion")
#'
#' # Alternatively, return a matrix
#' pairwise_ibs(example_gt, type = "proportion", as_matrix = TRUE)
#'
#' # Adjust block_size
#' pairwise_ibs(example_gt, block_size = 2)
#'
#' # Change type
#' pairwise_ibs(example_gt, type = "adjusted_counts")
#' pairwise_ibs(example_gt, type = "raw_counts")
#'
#' # PLINK 1.9's default correction for missing genotypes
#' pairwise_ibs(example_gt,
#'   type = "adjusted_counts",
#'   missing_correction = "allele_freq"
#' )
pairwise_ibs <- function(
  x,
  as_matrix = FALSE,
  type = c(
    "proportion",
    "adjusted_counts",
    "raw_counts"
  ),
  block_size = bigstatsr::block_size(nrow(x)),
  n_cores = bigstatsr::nb_cores(),
  missing_correction = c("flat", "allele_freq"),
  allele_freqs = NULL
) {
  stopifnot_diploid(x)
  type <- match.arg(type)
  X <- attr(x$genotypes, "fbm") # convenient pointer #nolint
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  ibs_matrix <- snp_ibs(
    X,
    ind.row = x_ind_row,
    ind.col = x_ind_col,
    type = type,
    block.size = block_size,
    n_cores = n_cores,
    missing_correction = missing_correction,
    allele_freqs = allele_freqs
  )
  if (inherits(ibs_matrix, "matrix")) {
    dimnames(ibs_matrix) <- list(x$id, x$id)
  } else {
    # else if we have a list of two count matrices
    attr(ibs_matrix[[1]], "indiv_names") <- x$id
    attr(ibs_matrix[[2]], "indiv_names") <- x$id
  }

  if (type == "proportion" || type == "adjusted_counts") {
    if (as_matrix) {
      # add a class to indicate this is a pairwise matrix
      class(ibs_matrix) <- c("pairwise_matrix", class(ibs_matrix))
      return(ibs_matrix)
    } else {
      return(tidy_dist_matrix(ibs_matrix))
    }
  }

  ibs_matrix
}
