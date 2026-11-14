#' Compute the Genetic Distances between individuals of a `gen_tibble` object
#'
#' This function computes the genetic distance between each pair of
#' individuals, as ngsDist (Vieira et al. 2016) does from genotypes: at each
#' locus genotyped in both individuals, a score of their two genotypes, 0 for
#' the same genotype, 0.5 for a homozygote and a heterozygote, and 1 for
#' opposite homozygotes, averaged over the loci, then corrected for
#' saturation (by default \eqn{-\log(1 - d)}{-log(1 - d)}, ngsDist's
#' default). The distances can be used to build trees, e.g. neighbour-joining
#' trees with the package `ape` (see the examples).
#'
#' Imputed genotypes are used if `x` is set to use them (see
#' [gt_set_imputed()]), and treated as missing otherwise.
#'
#' `pairwise_dist()` is an S3 generic. The method for a [`gen_tibble`]
#' computes the matrix with the method for its `genotypes` column (of class
#' `vctrs_bigSNP`), and labels it with the `id` of the individuals; grouping
#' is ignored. Other packages can add methods for their own classes of
#' genotypes (e.g. `tidypopgl`, from genotype likelihoods).
#'
#' @section The distance and its corrections:
#' With the default `score = "genotype"`, the uncorrected distance d of a
#' pair is the proportion of their alleles that are not shared identical by
#' state (at each locus, 0, 1 or 2 of their 2 alleles differ): \eqn{1 - d}{1 -
#' d} is the proportion of [pairwise_ibs()]. With `score = "nucleotide"`
#' (ngsDist's `--avg_nuc_dist`), two heterozygotes also score 0.5, and d is
#' the probability that an allele drawn at random from each individual
#' differs: \eqn{1 - d}{1 - d} is the allele sharing of
#' [pairwise_allele_sharing()] (whose diagonal, the allele sharing of an
#' individual with itself, is not 1: the distance of an individual to itself
#' is 0 here, as in ngsDist).
#'
#' A proportion is bounded (at most 1), so it saturates: the more two
#' individuals have diverged, the less each further change adds to it.
#' Tree-building methods such as neighbour joining and minimum evolution
#' assume that distances are additive (the distance between two individuals
#' is the sum of the lengths of the branches that join them), and a saturated
#' distance makes the long paths, and so the deep branches, too short. The
#' corrections stretch the large distances. All three have the form
#' \eqn{-b \log(1 - d / b)}{-b log(1 - d / b)}, where b is the ceiling that
#' d would reach between individuals that have diverged without limit:
#' * "log" (the default, ngsDist's): b = 1, \eqn{-\log(1 - d)}{-log(1 - d)},
#'   minus the logarithm of the proportion of alleles shared. The logarithm
#'   turns a similarity that falls by the same factor in each generation into
#'   a distance that grows linearly with time; it is the transformation of
#'   Reynolds et al.'s (1983) coancestry distance
#'   \eqn{-\log(1 - \theta)}{-log(1 - theta)}, which grows linearly with the
#'   generations since two populations split (divided by twice their size)
#'   when drift alone acts in populations of constant size. ngsDist's code
#'   gives this reason: distances additive under a constant population size,
#'   as neighbour joining and minimum evolution assume. For the alleles shared
#'   by two individuals it is a heuristic, not derived from a model of
#'   mutation: it only assumes that the alleles shared would fall to none
#'   with unlimited divergence.
#' * "jc": b = 3/4, Jukes and Cantor's (1969) correction for nucleotide
#'   sequences, where any base is equally likely to replace any other, so
#'   that unrelated sequences differ at 3/4 of their sites. It fits d as a
#'   divergence per nucleotide, over all the positions sequenced; on
#'   polymorphic loci alone, each with two alleles, its ceiling has no
#'   meaning, and d can reach 0.75 (`NA`, with a warning).
#' * "none": d itself (the p-distance; b infinite), the easiest to interpret.
#'
#' For small d the three agree: \eqn{-b \log(1 - d / b) \approx d + d^2 /
#' (2b)}{-b log(1 - d / b) ~ d + d^2 / (2b)}, so the log correction
#' increases a small d by a fraction of about d / 2, and a large d by more:
#' to 0.105 for d = 0.1 (5% more), 0.357 for d = 0.3 (19%), 0.693 for d = 0.5
#' (39%), and to infinity for d = 1 (opposite homozygotes at every locus;
#' `Inf`, with a warning). The corrections keep the order of the distances,
#' so heatmaps look alike and the nearest individuals stay the nearest; they
#' change the lengths of the branches of trees and, where distances are
#' large, can change their topology.
#'
#' @references Vieira, F. G., Lassalle, F., Korneliussen, T. S., Fumagalli, M.
#'   (2016). Improving the estimation of genetic distances from
#'   next-generation sequencing data. Biological Journal of the Linnean
#'   Society, 117(1), 139-149.
#'
#'   Reynolds, J., Weir, B. S., Cockerham, C. C. (1983). Estimation of the
#'   coancestry coefficient: basis for a short-term genetic distance.
#'   Genetics, 105(3), 767-779.
#'
#'   Jukes, T. H., Cantor, C. R. (1969). Evolution of protein molecules. In
#'   H. N. Munro (Ed.), Mammalian Protein Metabolism (pp. 21-132). Academic
#'   Press.
#' @param x a `gen_tibble` object, or a vector of class `vctrs_bigSNP` (the
#'   `genotypes` column of a [`gen_tibble`], whose individuals are then
#'   labelled with the names of the vector).
#' @param as_matrix boolean, determining whether the results should be a square
#'   symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)
#' @param score the scores of the pairs of genotypes: "genotype" (the
#'   default), 0 for the same genotype, 0.5 for a homozygote and a
#'   heterozygote, 1 for opposite homozygotes; or "nucleotide", the same but
#'   0.5 for two heterozygotes (ngsDist's `--avg_nuc_dist`).
#' @param correction the correction of the distances d for saturation
#'   (ngsDist's `--evol_model`; see "The distance and its corrections"):
#'   "log" (the default), \eqn{-\log(1 - d)}{-log(1 - d)}; "none", d itself
#'   (the p-distance); or "jc", Jukes and Cantor's \eqn{-3/4 \log(1 -
#'   4d/3)}{-3/4 log(1 - 4d/3)}, for d a divergence per nucleotide (`NA` for
#'   d of at least 0.75, with a warning).
#' @param block_size maximum number of loci read at once. Genotypes are
#'   compared in a compact bit-packed form, so this has little effect on
#'   memory use.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns a tidy tibble with one row per pair of individuals (columns
#'   `item1`, `item2`, `value`, the distance, and `n_loci`, the number of loci
#'   genotyped in both), or, with `as_matrix = TRUE`, a square matrix of the
#'   distances (0 on the diagonal) with the numbers of loci of each pair as
#'   attribute `n_loci`. Pairs without any locus genotyped in both have `NA`
#'   distances.
#' @seealso [pairwise_ibs()], [pairwise_allele_sharing()]
#' @rdname pairwise_dist
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # the distances, log corrected
#' pairwise_dist(example_gt)
#'
#' # the average nucleotide distance, uncorrected
#' pairwise_dist(example_gt, score = "nucleotide", correction = "none")
#'
#' # a neighbour-joining tree
#' dist_matrix <- pairwise_dist(example_gt, as_matrix = TRUE)
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   tree <- ape::nj(stats::as.dist(dist_matrix))
#'   plot(tree, type = "unrooted")
#' }
pairwise_dist <- function(x, ...) {
  UseMethod("pairwise_dist", x)
}

#' @export
#' @rdname pairwise_dist
pairwise_dist.tbl_df <- function(x, as_matrix = FALSE, ...) {
  stopifnot_has_genotypes(x)
  # hand the genotypes column over to its method (with its own arguments and
  # defaults), which returns a pairwise_matrix, and label it with the ids of
  # the tibble
  dist_matrix <- pairwise_dist(x$genotypes, as_matrix = TRUE, ...)
  dimnames(dist_matrix) <- list(x$id, x$id)
  if (as_matrix) {
    return(dist_matrix)
  }
  tidy_pairwise_dist(dist_matrix)
}

#' @export
#' @rdname pairwise_dist
pairwise_dist.vctrs_bigSNP <- function(
  x,
  as_matrix = FALSE,
  score = c("genotype", "nucleotide"),
  correction = c("log", "none", "jc"),
  block_size = bigstatsr::block_size(length(x)),
  n_cores = bigstatsr::nb_cores(),
  ...
) {
  rlang::check_dots_empty()
  stopifnot_diploid(x)
  score <- match.arg(score)
  correction <- match.arg(correction)
  X <- attr(x, "fbm") # nolint
  ind_row <- vctrs::vec_data(x)
  ind_col <- show_loci(x)$big_index
  if (length(ind_row) < 2) {
    stop("the distances need at least 2 individuals")
  }
  if (score == "genotype") {
    # the counts of alleles identical by state and of alleles compared (twice
    # the loci genotyped in both), which snp_ibs() returns as FBMs: made in a
    # private directory, deleted on exit
    local_private_fbm_dir()
    counts <- snp_ibs(
      X,
      ind.row = ind_row,
      ind.col = ind_col,
      type = "raw_counts",
      block.size = block_size,
      n_cores = n_cores
    )
    valid_n <- counts$valid_n[, , drop = FALSE]
    # the proportion of alleles that differ
    d <- 1 - counts$ibs[, , drop = FALSE] / valid_n
    n_loci <- valid_n / 2
  } else {
    # the sums of allele sharing: the probability that random alleles of the
    # two individuals differ is 1 - 0.5 * (1 + num / denom)
    counts <- .allele_sharing_counts(
      X, ind_row, ind_col, block_size, n_cores
    )
    d <- 0.5 * (1 - counts$num / counts$denom)
    n_loci <- counts$denom
  }
  # pairs without loci genotyped in both have no distance
  d[n_loci == 0] <- NA_real_
  d <- dist_correction(d, correction)
  diag(d) <- 0
  # the individuals are labelled with the names of the vector (the ids of the
  # gen_tibble)
  dimnames(d) <- list(names(x), names(x))
  attr(d, "n_loci") <- unname(n_loci)
  class(d) <- c("pairwise_matrix", class(d))
  if (as_matrix) {
    return(d)
  }
  tidy_pairwise_dist(d)
}

# The correction of the distances d for saturation (see pairwise_dist()), as
# ngsDist's gen_dist(); infinite (log) and undefined (jc) distances are
# given with a warning, where ngsDist gives them silently
dist_correction <- function(d, correction) {
  if (correction == "log") {
    if (any(d >= 1, na.rm = TRUE)) {
      warning("some distances are 1, so -log(1 - d) is infinite", call. = FALSE)
    }
    d <- -log(1 - d)
  } else if (correction == "jc") {
    beyond <- !is.na(d) & d >= 0.75
    if (any(beyond)) {
      warning(
        "some distances are at least 0.75: their Jukes-Cantor correction is NA",
        call. = FALSE
      )
    }
    d <- suppressWarnings(-log(1 - (d * 4 / 3)) * 3 / 4)
    d[beyond] <- NA_real_
  }
  d
}

# A tidy tibble of a matrix of distances (see tidy_dist_matrix()), with a
# column n_loci if the matrix has an attribute n_loci (the numbers of loci of
# each pair), so that methods of pairwise_dist() from other packages can give
# it too
tidy_pairwise_dist <- function(dist_matrix) {
  n_loci <- attr(dist_matrix, "n_loci")
  attr(dist_matrix, "n_loci") <- NULL
  class(dist_matrix) <- NULL
  dist_tbl <- tidy_dist_matrix(dist_matrix)
  if (!is.null(n_loci)) {
    # (the same pairs, in the same order)
    dimnames(n_loci) <- dimnames(dist_matrix)
    dist_tbl$n_loci <- tidy_dist_matrix(n_loci)$value
  }
  dist_tbl
}
