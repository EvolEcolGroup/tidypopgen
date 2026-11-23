#' Principal Component Analysis for `gen_tibble` objects
#'
#' Computes the leading principal components of the genotypes of a
#' `gen_tibble`, with one of two methods. `method = "partial"` (the default)
#' is exact and handles missing genotypes as `plink --pca` (PLINK 1.9) does,
#' without imputation; `method = "random"` needs complete or imputed
#' genotypes, but much less memory, for very large numbers of individuals.
#' For datasets with linkage disequilibrium, [gt_pca_autoSVD()] implements the
#' algorithm of `bigsnpr`, which removes regions of long-range LD and outliers.
#'
#' @section Methods: Both methods compute the leading eigenvectors of the
#'   genetic relationship matrix (GRM) of the standardised genotypes
#'   \eqn{z_{il}} (by default \eqn{z_{il} = (g_{il} - 2p_l) /
#'   \sqrt{2p_l(1-p_l)}}, with \eqn{g_{il}} the number of alternate alleles
#'   and \eqn{p_l} their frequency):
#'
#'   * `"partial"` computes the GRM, then its leading eigenvectors (an exact,
#'     truncated decomposition). It is the method of choice up to some tens of
#'     thousands of individuals: the GRM takes 8 bytes per pair of
#'     individuals (0.2 GB for 5,000 individuals, 3.2 GB for 20,000), and the
#'     time grows with the square of the number of individuals and linearly
#'     with the number of loci.
#'   * `"random"` runs [bigstatsr::big_randomSVD()], which never builds the
#'     GRM: it multiplies blocks of the genotypes by vectors until the
#'     components converge (to `tol`), with memory linear in the numbers of
#'     individuals and loci. It is slower than `"partial"` (e.g. 175 s
#'     against 31 s on 1,043 individuals and 643,733 loci), but works with
#'     any number of individuals.
#'
#'   With complete (or imputed) genotypes, both give the singular value
#'   decomposition of the standardised genotypes, and `"partial"` gives the
#'   same results as the former [gt_pca_partialSVD()].
#'
#' @section Missing genotypes: With `method = "partial"`, missing genotypes
#'   do not need to be imputed: as in PLINK 1.9 (`plink --make-rel`, which
#'   `plink --pca` decomposes), each pair of individuals is compared over the
#'   loci genotyped in both,
#'   \deqn{G_{ij} = \frac{1}{n_{ij}} \sum_l z_{il} z_{jl},}
#'   with \eqn{z_{il} = 0} when \eqn{g_{il}} is missing, \eqn{n_{ij}} the
#'   number of loci genotyped in both individuals, and allele frequencies
#'   from the genotyped individuals. This keeps individuals with many missing
#'   genotypes (e.g. ancient DNA) in place, whereas imputation pulls them
#'   towards the centre of the PCA, or makes missingness an axis of its own.
#'   The results match those of PLINK 1.9, also with very high missingness.
#'   Imputed genotypes are used only if the `gen_tibble` is set to use them
#'   (with [gt_set_imputed()]). Pairs of individuals with no locus genotyped
#'   in both have no defined relationship and give an error: remove the
#'   individuals with the most missing genotypes (see [indiv_missingness()]).
#'   With missing genotypes, the GRM is not necessarily positive
#'   semi-definite; only its largest (positive) eigenvalues are used.
#'
#'   With `method = "random"`, the genotypes must be complete or imputed (see
#'   [gt_impute_simple()]); imputed genotypes are used automatically when
#'   they exist.
#'
#'   The PCA records which genotypes it used (see Value). With missing
#'   genotypes, [gt_pcadapt()] tests each locus on the individuals genotyped
#'   at it, and [gt_dapc()] estimates the loadings of the loci without
#'   counting missing genotypes as 0 (for which it needs the `gen_tibble`).
#'
#' @section Scaling: `fun_scaling = NULL` (the default) uses the binomial
#'   scaling above, as [bigsnpr::snp_scaleBinom()] (and PLINK) do, with the
#'   allele frequencies of the genotyped individuals. A function that returns
#'   a centre and a scale for each locus, such as
#'   [bigsnpr::snp_scaleAlpha()] or [bigstatsr::big_scale()], can be used
#'   instead. A scaling computed from the genotypes (such as these) needs
#'   complete or imputed genotypes, and so does any scaling with `"random"`
#'   (one that does not depend on the genotypes, e.g.
#'   `big_scale(center = FALSE, scale = FALSE)`, works with missing genotypes
#'   with `"partial"`). With `"partial"`, loci
#'   with a zero scale (monomorphic loci) add nothing to the GRM, but are
#'   counted in \eqn{n_{ij}}, as in PLINK (filter them out first); with
#'   `"random"`, they give an error.
#'
#' @section Pseudohaploids: Pseudohaploid individuals (see
#'   [gt_pseudohaploid()]) are coded 0 or 2, and enter the PCA as homozygous
#'   diploids, as in PLINK and `smartpca` with pseudohaploid ancient DNA:
#'   half their dosage is their allele frequency, as for diploids, while the
#'   allele frequencies of the default scaling count one allele per
#'   pseudohaploid (as [loci_alt_freq()] does). Their own entries on the
#'   diagonal of the GRM are about twice those of diploids, which does not
#'   change the relationships between individuals.
#'
#' @section Speed: With `method = "partial"`, the matrix products use either
#'   the BLAS library that R uses (`engine = "blas"`), which is fastest when
#'   it is an optimised one (OpenBLAS, MKL, Apple Accelerate), or the
#'   vectorised products of the Eigen library in parallel on `n_cores`
#'   (`engine = "internal"`), which do not depend on the BLAS and are many
#'   times faster than R's reference BLAS (the default on Windows and in the
#'   CRAN builds of R for macOS). The internal products are in single
#'   precision for blocks of up to 1,024 loci, summed in double precision
#'   (twice as fast; the relative error of the GRM is about 1e-7, well below
#'   the 6 significant digits that PLINK writes). `engine = "auto"` times both
#'   on a small product, once per session, and uses the faster. The number of
#'   loci genotyped in both individuals is counted on bit-packed genotypes, as
#'   in PLINK. On a laptop with 4 cores and OpenBLAS, 4,000 individuals and
#'   100,000 loci with 45% missing genotypes took 32 s (62 s with PLINK 1.9),
#'   and 1,000 individuals and 300,000 complete loci 7 s (8 s with PLINK, 11 s
#'   with the former [gt_pca_partialSVD()]); with R's reference BLAS, 39 s and
#'   12 s.
#'
#' @param x a `gen_tbl` object, with diploid (or pseudohaploid) genotypes
#' @param k the number of principal components to compute. Default is `10`.
#' @param method `"partial"` (the default) or `"random"` (see Methods).
#' @param fun_scaling `NULL` (the default) for the binomial scaling, or a
#'   scaling function (see Scaling).
#' @param n_cores the number of cores. With `"partial"`, they read the
#'   genotypes, count the loci genotyped in both individuals and run the
#'   internal matrix products; the BLAS uses its own threads (see the
#'   `RhpcBLASctl` package to set them). With `"random"`, they run the
#'   matrix products of [bigstatsr::big_prodVec()] and
#'   [bigstatsr::big_cprodVec()] (in one process); with other `fun_prod` or
#'   `fun_cprod`, the parallel SVD of [bigstatsr::big_randomSVD()] runs on
#'   `n_cores - 1` workers and one more process that coordinates them (or on
#'   one core, with fewer than 3 cores).
#' @param block_size for `"partial"`, the number of loci decoded at a time
#'   (each block takes 8 bytes per individual and locus). By default, blocks
#'   of about 128 MB (between 1,024 and 16,384 loci).
#' @param engine for `"partial"`, the matrix products: `"auto"` (the
#'   default), `"blas"` or `"internal"` (see Speed).
#' @param tol for `"random"`, the precision parameter of [RSpectra::svds()].
#'   Default is `1e-4`.
#' @param verbose for `"random"`, should some progress be printed? Default is
#'   `FALSE`.
#' @param fun_prod for `"random"`, a function that takes 6 arguments (in this
#'   order): a matrix-like object `X`, a vector `x`, a vector of row indices
#'   `ind.row` of `X`, a vector of column indices `ind.col` of `X`, a vector
#'   of column centers (corresponding to `ind.col`), and a vector of column
#'   scales (corresponding to `ind.col`), and computes the product of `X`
#'   (subsetted and scaled) with `x`. Default is [bigstatsr::big_prodVec()].
#' @param fun_cprod for `"random"`, the same as `fun_prod`, but for the
#'   *transpose* of `X`. Default is [bigstatsr::big_cprodVec()].
#' @param total_var a boolean indicating whether to compute the total
#'   variance, needed for the proportion of variance explained by each
#'   component (see [`gt_pca_tidiers`]). Default is `TRUE`.
#' @returns a `gt_pca` object, which is a subclass of `big_SVD`; this is an S3
#'   list with elements:
#' - `d`, the singular values (the eigenvalues of the GRM are `d^2 / m`, with
#'   `m` the number of loci, as given by `plink --pca`)
#' - `u`, the scores for each sample on each component (the left singular
#'   vectors, i.e. the eigenvectors of the GRM)
#' - `v`, the loadings (the right singular vectors; with missing genotypes
#'   and `"partial"`, \eqn{Z^T u / d}, with missing genotypes as 0)
#' - `center`, the centering vector
#' - `scale`, the scaling vector
#' - `method`, `"partial"` or `"random"`
#' - `call`, the call to the function
#' - `loci`, the loci used
#' - `square_frobenius`, the total variance (with `total_var = TRUE`; with
#'   missing genotypes and `"partial"`, the number of loci times the trace of
#'   the GRM)
#'
#' Two attributes record the genotypes that the PCA was computed on:
#' - `imputed`, `TRUE` if the PCA used imputed genotypes
#' - `missing`, `TRUE` if the genotypes had missing values (only possible
#'   with `"partial"`)
#'
#' PCAs made by earlier versions of `tidypopgen` do not have these
#' attributes; their methods needed complete or imputed genotypes.
#'
#' @note Rather than accessing elements of a `gt_pca` object directly, it is
#'   better to use `tidy` and `augment`. See [`gt_pca_tidiers`].
#' @references Purcell S. et al. (2007) PLINK: a tool set for whole-genome
#'   association and population-based linkage analyses. American Journal of
#'   Human Genetics 81: 559-575.
#'
#'   Chang C.C. et al. (2015) Second-generation PLINK: rising to the challenge
#'   of larger and richer datasets. GigaScience 4: 7.
#'
#'   Privé F., Luu K., Blum M.G.B., McGrath J.J. and Vilhjálmsson B.J. (2020)
#'   Efficient toolkit implementing best practices for principal component
#'   analysis of population genetic data. Bioinformatics 36: 4449-4457.
#' @seealso [pairwise_grm()] with `method = "plink"` for the GRM itself;
#'   [gt_pca_autoSVD()] for data with linkage disequilibrium.
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' # the example has missing genotypes, which are compared pairwise
#' gt_pca(example_gt, k = 2)
#'
#' # three simulated populations, with up to 90% missing genotypes
#' set.seed(1)
#' sim_gt <- gt_simulate_bn(
#'   n_indiv = c(20, 20, 20), n_loci = 2000, fst = 0.1,
#'   missing = rep(seq(0, 0.9, length.out = 20), 3)
#' )
#' sim_pca <- gt_pca(sim_gt, k = 2, n_cores = 1)
#' autoplot(sim_pca, type = "scores")
#'
#' # the random method needs complete or imputed genotypes, without
#' # monomorphic loci
#' sim_gt <- select_loci_if(sim_gt, loci_maf(genotypes) > 0)
#' sim_gt <- gt_impute_simple(sim_gt, method = "mode")
#' gt_pca(sim_gt, k = 2, method = "random")
#' @family gt_pca_functions
gt_pca <- function(x, k = 10, method = c("partial", "random"),
                   fun_scaling = NULL, n_cores = bigstatsr::nb_cores(),
                   block_size = NULL, engine = c("auto", "blas", "internal"),
                   tol = 1e-4, verbose = FALSE,
                   fun_prod = bigstatsr::big_prodVec,
                   fun_cprod = bigstatsr::big_cprodVec,
                   total_var = TRUE) {
  stopifnot_dip_pseudo(x)
  method <- match.arg(method)
  engine <- match.arg(engine)
  n_ind <- nrow(x)
  if (k < 1 || k > n_ind - 1) {
    stop("'k' should be between 1 and the number of individuals - 1")
  }
  if (method == "random") {
    # (pca_random_svd() uses imputed genotypes if they exist; otherwise they
    # must be complete, which its scaling checks)
    # (bigsnpr's binomial scaling counts two alleles per individual, so
    # pseudohaploids get the frequencies of loci_alt_freq(), as with
    # "partial"; that scaling ignores missing genotypes, so all of them are
    # then looked for, as for a custom scaling)
    # (computed after switching on the imputed genotypes, if any, as
    # pca_random_svd() does before the default scaling, so that automatic and
    # explicit activation give the same frequencies)
    pseudo_scaling <- is.null(fun_scaling) && is_pseudohaploid(x)
    if (pseudo_scaling) {
      if (gt_has_imputed(x) && gt_uses_imputed(x) == FALSE) {
        gt_set_imputed(x, set = TRUE)
        on.exit(gt_set_imputed(x, set = FALSE), add = TRUE)
      }
      fun_scaling <- pseudohaploid_scaling(x, n_cores)
    }
    this_pca <- pca_random_svd(x,
      k = k,
      fun_scaling = if (is.null(fun_scaling)) {
        bigsnpr::snp_scaleBinom()
      } else {
        fun_scaling
      },
      tol = tol, verbose = verbose, n_cores = n_cores, fun_prod = fun_prod,
      fun_cprod = fun_cprod, total_var = total_var,
      what = "gt_pca(method = \"random\")",
      hint = "or use method = \"partial\"",
      custom_scaling = !is.null(fun_scaling)
    )
    this_pca$method <- "random"
    this_pca$call <- match.call()
    return(this_pca)
  }

  # method "partial": the GRM, then its leading eigenvectors
  # (imputed genotypes are used only if x is set to use them)
  imputed <- gt_has_imputed(x) && gt_uses_imputed(x)
  if (is.null(block_size)) block_size <- grm_block_size(n_ind)
  use_blas <- grm_use_blas(engine, n_cores)
  scaling <- pca_scaling(x, fun_scaling, n_cores)
  col_ind <- show_loci(x)$big_index
  m <- length(col_ind)
  grm <- grm_plink(x,
    center = scaling$center, scale = scaling$scale, n_cores = n_cores,
    block_size = block_size, use_blas = use_blas
  )
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
    BM = .gt_get_fbm(x), rowInd = .gt_fbm_rows(x), colInd = col_ind,
    center = scaling$center, scale = scaling$scale, U = u, d = d,
    block_size = block_size, ncores = n_cores, use_blas = use_blas
  )
  v <- matrix_int_names(v)
  row_names(v) <- loci_names(x)
  this_pca <- list(
    d = d, u = u, v = v, center = scaling$center, scale = scaling$scale,
    method = "partial", call = match.call(), loci = show_loci(x)
  )
  if (total_var) this_pca$square_frobenius <- m * sum(diag(grm))
  class(this_pca) <- c("gt_pca", "big_SVD")
  # record the genotypes the PCA was computed on: imputed or not, and with
  # missing genotypes or not (counted by grm_plink_cpp()), for the functions
  # that treat the two differently (gt_dapc() and gt_pcadapt())
  pca_set_genotypes(this_pca,
    imputed = imputed,
    missing = attr(grm, "n_missing") > 0
  )
}

# Record, as attributes of a gt_pca, whether the PCA used imputed genotypes
# ("imputed") and whether the genotypes it decomposed had missing values
# ("missing"; only possible with method "partial")
pca_set_genotypes <- function(pca, imputed, missing) {
  attr(pca, "imputed") <- imputed
  attr(pca, "missing") <- missing
  pca
}

# Centering and scaling of the loci for method "partial": the binomial
# scaling with the allele frequencies of the genotyped individuals
# (fun_scaling = NULL), or those of a scaling function (which needs complete
# or imputed genotypes). Loci with no genotype, or no scale, get a scale of
# 0, so that they add nothing to the GRM.
pca_scaling <- function(x, fun_scaling, n_cores) {
  if (is.null(fun_scaling)) {
    # on the genotypes column, so that a grouped gen_tibble gives the
    # frequencies of all individuals (not one column per group)
    freq <- loci_alt_freq(x$genotypes, n_cores = n_cores)
    center <- 2 * freq
    scale <- sqrt(2 * freq * (1 - freq))
  } else {
    # bigstatsr's check for missing values samples the whole backing matrix
    # (including individuals and loci that x does not use): turn it off;
    # missing genotypes give missing values in the scaling, checked below
    local_no_bigstatsr_checks()
    scaling <- tryCatch(
      fun_scaling(.gt_get_fbm(x),
        ind.row = .gt_fbm_rows(x),
        ind.col = .gt_fbm_cols(x)
      ),
      error = function(e) {
        # (bigstatsr's own check for missing values is off, so scalings can
        # fail on them with a cryptic error: count them, only now)
        if (any_missing_genotypes(x, n_cores = n_cores)) {
          stop_missing_genotypes("gt_pca() with a custom 'fun_scaling'",
            hint = "or use the default scaling (fun_scaling = NULL)",
            imputed_auto = FALSE
          )
        }
        stop("'fun_scaling' failed: ", conditionMessage(e), call. = FALSE)
      }
    )
    center <- scaling$center
    scale <- scaling$scale
    # missing genotypes give missing values in the scaling (an exact check)
    if (anyNA(center) || anyNA(scale)) {
      stop_missing_genotypes("gt_pca() with a custom 'fun_scaling'",
        hint = "or use the default scaling (fun_scaling = NULL)",
        imputed_auto = FALSE
      )
    }
  }
  center[!is.finite(center)] <- 0
  scale[!is.finite(scale)] <- 0
  list(center = center, scale = scale)
}

# The GRM of method "partial" (and of pairwise_grm(method = "plink")): the
# standardised genotypes are compared over the loci genotyped in both
# individuals, as in PLINK 1.9 (--make-rel). Returns an n x n matrix, NaN for
# pairs with no locus genotyped in both.
grm_plink <- function(x, center, scale, n_cores, block_size, use_blas) {
  code256 <- .gt_get_fbm(x)$code256
  # only genotypes 0, 1 and 2 (and missing) are supported, not dosages
  if (any(is.finite(code256) & !code256 %in% c(0, 1, 2))) {
    stop(
      "the PLINK relationship matrix needs genotypes 0, 1 and 2, ",
      "not dosages"
    )
  }
  grm_plink_cpp(
    BM = .gt_get_fbm(x), rowInd = .gt_fbm_rows(x),
    colInd = show_loci(x)$big_index, center = center, scale = scale,
    block_size = block_size, ncores = n_cores, use_blas = use_blas
  )
}

# The binomial scaling of the loci of a gen_tibble with pseudohaploids, for
# method "random" (bigsnpr::snp_scaleBinom() would count two alleles per
# pseudohaploid): centre 2 f and scale sqrt(2 f (1 - f)), with f from
# loci_alt_freq(), as a scaling function of bigstatsr (the values of the
# columns ind.col of the backing matrix). Loci without a scale get 1, so
# that the SVD does not divide by 0 (their centred genotypes are 0).
pseudohaploid_scaling <- function(x, n_cores) {
  freq <- loci_alt_freq(x$genotypes, n_cores = n_cores)
  fbm_cols <- show_loci(x)$big_index
  scale <- sqrt(2 * freq * (1 - freq))
  scale[!is.finite(scale) | scale == 0] <- 1
  function(X, ind.row, ind.col) { # nolint: object_name_linter
    pos <- match(ind.col, fbm_cols)
    data.frame(center = 2 * freq[pos], scale = scale[pos])
  }
}

# Method "random": bigstatsr::big_randomSVD(), using imputed genotypes if they
# exist (also used by the deprecated gt_pca_randomSVD()); `what` and `hint`
# are for the error on missing genotypes (see stop_missing_genotypes()), and
# `custom_scaling` says whether the user gave fun_scaling (see
# scaling_lookup())
pca_random_svd <- function(x, k, fun_scaling, tol, verbose, n_cores,
                           fun_prod, fun_cprod, total_var,
                           what = "gt_pca_randomSVD()", hint = NULL,
                           custom_scaling = TRUE) {
  if (gt_has_imputed(x) && gt_uses_imputed(x) == FALSE) {
    gt_set_imputed(x, set = TRUE)
    on.exit(gt_set_imputed(x, set = FALSE), add = TRUE)
  }
  if (n_cores > 1) {
    # Remove checking for two levels of parallelism
    .old_opt <- getOption("bigstatsr.check.parallel.blas", TRUE)
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = .old_opt), add = TRUE)
  }
  # the SVD needs complete genotypes: bigstatsr's check samples the whole
  # backing matrix (including individuals and loci that x does not use), so
  # turn it off; the scaling, computed once for all the loci, checks the
  # genotypes of x (see scaling_lookup())
  # (keeping bigstatsr's check of n_cores)
  local_no_bigstatsr_checks(n_cores = n_cores)
  fun_scaling <- scaling_lookup(x, fun_scaling,
    what = what, hint = hint, n_cores = n_cores,
    custom_scaling = custom_scaling
  )

  X <- attr(x$genotypes, "fbm") # convenient pointer #nolint
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)

  # use at most n_cores cores: the parallel SVD of big_randomSVD() runs
  # 1 + ncores processes (one runs RSpectra, ncores compute the products).
  # With bigstatsr's own products (the default), run its sequential SVD,
  # with the products on n_cores threads in one process (on HGDP650, 1,043
  # individuals and 488,797 loci, 4 cores: 15.1 s, against 17.2 s for the
  # parallel SVD on 1 + 4 processes and 20.9 s on 1 + 3); with products
  # given by the user, run the parallel SVD on n_cores - 1 workers (n_cores
  # processes in all), or the sequential one with fewer than 3 cores
  svd_cores <- svd_n_cores(n_cores)
  default_prods <- identical(fun_prod, bigstatsr::big_prodVec) &&
    identical(fun_cprod, bigstatsr::big_cprodVec)
  if (default_prods) {
    svd_cores <- 1
    prod_cores <- n_cores
    fun_prod <- function(X, x, ind.row, ind.col, center, scale) { # nolint
      bigstatsr::big_prodVec(X, x, ind.row, ind.col,
        center = center, scale = scale, ncores = prod_cores
      )
    }
    fun_cprod <- function(X, x, ind.row, ind.col, center, scale) { # nolint
      bigstatsr::big_cprodVec(X, x, ind.row, ind.col,
        center = center, scale = scale, ncores = prod_cores
      )
    }
  }

  # big_randomSVD() makes a temporary FBM when run in parallel: keep it in a
  # private directory, deleted on exit (see with_private_fbm_dir())
  this_svd <- with_private_fbm_dir(bigstatsr::big_randomSVD(
    X,
    k = k,
    ind.row = .gt_fbm_rows(x),
    ind.col = .gt_fbm_cols(x),
    # (which also turns off bigstatsr's checks in the workers)
    fun.scaling = fun_scaling,
    tol = tol,
    verbose = verbose,
    ncores = svd_cores,
    fun.prod = fun_prod,
    fun.cprod = fun_cprod
  ))
  # add names to the scores (to match them to data later)
  row.names(this_svd$u) <- x$id
  this_svd$v <- matrix_int_names(this_svd$v)
  row_names(this_svd$v) <- loci_names(x)
  this_svd$method <- "randomSVD"
  this_svd$loci <- show_loci(x)
  class(this_svd) <- c("gt_pca", class(this_svd))
  if (total_var) {
    this_svd$square_frobenius <- square_frobenius(
      X,
      x_ind_row,
      x_ind_col,
      center = this_svd$center,
      scale = this_svd$scale
    )
  }
  # the genotypes were complete (or imputed, if x has imputed genotypes)
  pca_set_genotypes(this_svd, imputed = gt_has_imputed(x), missing = FALSE)
}

# The number of workers for the parallel randomised SVD of bigstatsr
# (big_randomSVD(), also in bigsnpr's snp_autoSVD()) that uses at most n_cores
# cores: it runs 1 + ncores processes (one runs RSpectra, ncores compute the
# products), so n_cores - 1 workers; with fewer than 3 cores, 1 (its
# sequential SVD, in one process)
# @param n_cores the number of cores asked for
# @returns the `ncores` to give to bigstatsr
# @keywords internal
# @noRd
svd_n_cores <- function(n_cores) {
  if (n_cores >= 3) n_cores - 1 else 1
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
