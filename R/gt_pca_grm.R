#' PCA for `gen_tibble` objects with missing genotypes, as in PLINK 1.9
#'
#' This function performs Principal Component Analysis on a `gen_tibble` with
#' missing genotypes, without imputation, by the eigen decomposition of a
#' genetic relationship matrix (GRM) that accounts for missing data as
#' `plink --pca` (PLINK 1.9) does: each pair of individuals is compared over
#' the loci genotyped in both. The results match those of PLINK 1.9, also with
#' very high missingness, as in ancient DNA.
#'
#' The GRM is the variance-standardised relationship matrix of PLINK 1.9
#' (`--make-rel`):
#' \deqn{G_{ij} = \frac{1}{n_{ij}} \sum_l z_{il} z_{jl}, \qquad
#' z_{il} = \frac{g_{il} - 2p_l}{\sqrt{2p_l(1-p_l)}},}
#' where \eqn{g_{il}} is the number of alternate alleles, \eqn{p_l} the
#' alternate allele frequency among the genotyped individuals, \eqn{z_{il} = 0}
#' when \eqn{g_{il}} is missing, and \eqn{n_{ij}} is the number of loci
#' genotyped in both individuals. Monomorphic loci add nothing to the sum but
#' are counted in \eqn{n_{ij}}, as in PLINK (which warns about them: filter
#' loci by minor allele frequency first). Dividing each pair by the loci it
#' shares keeps individuals with many missing genotypes in place; imputation
#' instead pulls them towards the centre of the PCA, or makes missingness an
#' axis of its own.
#'
#' This is the matrix returned by `pairwise_grm(x, method = "plink")`. It is
#' not the default of [pairwise_grm()], the estimator of Weir and Goudet
#' (2017), which is based on allele sharing: it gives the same weight to every
#' locus and measures relatedness relative to the average pair, whereas the
#' variance-standardised matrix weights loci by \eqn{1 / (2p(1-p))}, so that
#' rare variants count more, and is the matrix whose eigenvectors are the
#' principal components of PLINK and of the other `gt_pca` functions.
#'
#' The principal components are the leading eigenvectors of \eqn{G} (`u`),
#' and its eigenvalues are `d^2 / m`, with `m` the number of loci, as for the
#' other `gt_pca` functions (the eigenvalues of `plink --pca`). With missing
#' genotypes, \eqn{G} is not necessarily positive semi-definite; only its
#' largest (positive) eigenvalues are used. The loadings `v` are
#' \eqn{Z^T u / d}, with missing genotypes as 0, so that for complete data
#' they are the right singular vectors of the standardised genotypes, and the
#' result matches that of [gt_pca_partialSVD()].
#'
#' Pairs of individuals with no locus genotyped in both have no defined
#' relationship, and give an error: remove the individuals with the most
#' missing genotypes (see [indiv_missingness()]).
#'
#' The numerator is accumulated with the BLAS over blocks of `block_size`
#' loci, and the number of loci genotyped in both individuals is counted on
#' bit-packed genotypes, as in PLINK. The time grows with the square of the
#' number of individuals and linearly with the number of loci; the GRM is
#' kept in memory (8 bytes per pair of individuals, about 0.5 GB for 8,000
#' individuals). The products can use the BLAS library or Eigen (see
#' `engine`). Imputed genotypes are used if the `gen_tibble` is set to use
#' them (see [gt_set_imputed()]); otherwise they are treated as missing.
#'
#' @param x a `gen_tbl` object, with diploid genotypes
#' @param k the number of principal components to compute. Default is `10`.
#' @param n_cores the number of cores used to read the genotypes and count the
#'   loci genotyped in both individuals. The matrix products use the threads
#'   of the BLAS library (see the `RhpcBLASctl` package to set them).
#' @param block_size the number of loci decoded at a time (each block takes
#'   8 bytes per individual and locus). By default, blocks of about 128 MB
#'   (between 1,024 and 16,384 loci): larger blocks make the matrix products
#'   more efficient.
#' @param engine the matrix products: `"blas"` uses the BLAS library that R
#'   uses, which is fastest when it is an optimised one (OpenBLAS, MKL, Apple
#'   Accelerate); `"internal"` uses the vectorised products of the Eigen
#'   library in parallel on `n_cores`, which do not depend on the BLAS and are
#'   many times faster than R's reference BLAS (the default on Windows and in
#'   the CRAN builds of R for macOS). The internal products are in single
#'   precision for blocks of up to 1,024 loci, summed in double precision
#'   (twice as fast; the relative error of the GRM is about 1e-7, well below
#'   the 6 significant digits that PLINK writes). `"auto"` (the default) times
#'   both on a small product, once per session, and uses the faster.
#' @returns a `gt_pca` object, which is a subclass of `big_SVD`; this is an S3
#'   list with elements:
#' - `d`, the singular values (the eigenvalues of the GRM are `d^2 / m`)
#' - `u`, the scores for each sample on each component (the eigenvectors of
#'   the GRM)
#' - `v`, the loadings
#' - `center`, the centering vector (twice the alternate allele frequencies)
#' - `scale`, the scaling vector (\eqn{\sqrt{2p(1-p)}})
#' - `method`, `"grm"`
#' - `call`, the call to the function
#' - `loci`, the loci used
#' - `square_frobenius`, `m` times the trace of the GRM, so that the
#'   proportion of variance explained by a component is its eigenvalue
#'   divided by the trace of the GRM
#' @seealso [pairwise_grm()] with `method = "plink"` for the relationship
#'   matrix itself; [gt_pca_partialSVD()] for complete or imputed data.
#' @references Purcell S. et al. (2007) PLINK: a tool set for whole-genome
#'   association and population-based linkage analyses. American Journal of
#'   Human Genetics 81: 559-575.
#'
#'   Chang C.C. et al. (2015) Second-generation PLINK: rising to the challenge
#'   of larger and richer datasets. GigaScience 4: 7.
#'
#'   Weir B.S. and Goudet J. (2017) A unified characterization of population
#'   structure and relatedness. Genetics 206: 2085-2103.
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' # the example has missing genotypes, which gt_pca_partialSVD() can not use
#' gt_pca_grm(example_gt, k = 2)
gt_pca_grm <- function(x, k = 10, n_cores = bigstatsr::nb_cores(),
                       block_size = NULL,
                       engine = c("auto", "blas", "internal")) {
  stopifnot_diploid(x)
  if (is.null(block_size)) block_size <- grm_block_size(nrow(x))
  use_blas <- grm_use_blas(match.arg(engine), n_cores)
  geno_fbm <- .gt_get_fbm(x)
  n_ind <- nrow(x)
  if (k < 1 || k > n_ind - 1) {
    stop("'k' should be between 1 and the number of individuals - 1")
  }
  col_ind <- show_loci(x)$big_index
  m <- length(col_ind)
  grm <- grm_plink(x, n_cores = n_cores, block_size = block_size,
    use_blas = use_blas
  )
  center <- attr(grm, "center")
  scale <- attr(grm, "scale")
  # pairs without a locus genotyped in both have no relationship
  n_empty <- sum(is.nan(grm[lower.tri(grm)]))
  if (n_empty > 0) {
    stop(
      n_empty, " pairs of individuals have no locus genotyped in both; ",
      "remove the individuals with the most missing genotypes"
    )
  }
  # the largest eigenvalues (with missing genotypes, some can be negative)
  if (k > n_ind / 4) {
    eig <- eigen(grm, symmetric = TRUE)
    eig$values <- eig$values[seq_len(k)]
    eig$vectors <- eig$vectors[, seq_len(k), drop = FALSE]
  } else {
    eig <- RSpectra::eigs_sym(grm, k, which = "LA")
  }
  if (any(eig$values <= 0)) {
    stop(
      "only ", sum(eig$values > 0), " eigenvalues of the GRM are positive; ",
      "use a smaller 'k'"
    )
  }
  d <- sqrt(eig$values * m)
  u <- eig$vectors
  row.names(u) <- x$id
  # loadings: Z^T u / d, with missing genotypes as 0
  v <- grm_loadings_cpp(
    BM = geno_fbm, rowInd = .gt_fbm_rows(x), colInd = col_ind,
    center = center, scale = scale, U = u, d = d,
    block_size = block_size, ncores = n_cores, use_blas = use_blas
  )
  v <- matrix_int_names(v)
  row_names(v) <- loci_names(x)
  this_pca <- list(
    d = d, u = u, v = v, center = center, scale = scale,
    method = "grm", call = match.call(), loci = show_loci(x),
    square_frobenius = m * sum(diag(grm))
  )
  class(this_pca) <- c("gt_pca", "big_SVD")
  this_pca
}

# The variance-standardised relationship matrix of PLINK 1.9 (--make-rel),
# with missing genotypes compared over the loci genotyped in both individuals
# (see gt_pca_grm()); used by gt_pca_grm() and pairwise_grm(method = "plink").
# Returns an n x n matrix (NaN for pairs with no locus genotyped in both),
# with the centering (2p) and scaling (sqrt(2p(1 - p))) vectors as attributes.
grm_plink <- function(x, n_cores, block_size, use_blas) {
  geno_fbm <- .gt_get_fbm(x)
  code256 <- geno_fbm$code256
  # only genotypes 0, 1 and 2 (and missing) are supported, not dosages
  if (any(is.finite(code256) & !code256 %in% c(0, 1, 2))) {
    stop(
      "the PLINK relationship matrix needs genotypes 0, 1 and 2, ",
      "not dosages"
    )
  }
  # alternate allele frequencies among the genotyped individuals
  freq <- loci_alt_freq(x, n_cores = n_cores)
  center <- 2 * freq
  scale <- sqrt(2 * freq * (1 - freq))
  # loci with no genotype have no frequency: they add nothing
  center[!is.finite(center)] <- 0
  scale[!is.finite(scale)] <- 0
  grm <- grm_plink_cpp(
    BM = geno_fbm, rowInd = .gt_fbm_rows(x),
    colInd = show_loci(x)$big_index, center = center, scale = scale,
    block_size = block_size, ncores = n_cores, use_blas = use_blas
  )
  attr(grm, "center") <- center
  attr(grm, "scale") <- scale
  grm
}

# Whether the products of the PLINK relationship matrix use the BLAS (TRUE)
# or Eigen (FALSE). For engine "auto", both are timed on a small product
# (512 individuals, 512 loci) and the faster is used; the choice is kept for
# the session, for each number of cores.
grm_use_blas <- function(engine, n_cores) {
  if (engine != "auto") {
    return(engine == "blas")
  }
  key <- paste0("grm_use_blas_", n_cores)
  if (is.null(tidypopgen_cache[[key]])) {
    # the best of two runs of each, to discount the first call's overheads
    secs <- sapply(c(TRUE, FALSE), function(use_blas) {
      min(replicate(2, grm_rank_update_seconds(512L, 512L, use_blas, n_cores)))
    })
    tidypopgen_cache[[key]] <- secs[1] <= secs[2]
  }
  tidypopgen_cache[[key]]
}

# the session cache of grm_use_blas()
tidypopgen_cache <- new.env(parent = emptyenv())

# Default number of loci per block of the PLINK relationship matrix: about
# 128 MB of doubles, between 1,024 and 16,384 loci
grm_block_size <- function(n_ind) {
  as.integer(min(16384, max(1024, floor(2^27 / (8 * max(n_ind, 1))))))
}
