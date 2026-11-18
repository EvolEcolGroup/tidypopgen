#' Predict scores of a PCA
#'
#' Predict the PCA scores for a [`gt_pca`], either for the original data or
#' projecting new data.
#'
#' New individuals are projected on the loadings of the PCA, with their
#' genotypes centred and scaled as those of the PCA. The methods differ in how
#' they deal with missing genotypes:
#'
#' * `"least_squares"` (the default): for each individual, the least squares
#'   fit of its standardised genotypes on the loadings of the loci genotyped
#'   in that individual (as `lsqproject` in smartpca). Individuals with many
#'   missing genotypes are not shrunk towards the centre of the PCA. With
#'   complete genotypes, it gives the same scores as `"simple"`. Individuals
#'   with too few genotyped loci to fit the components get `NA`, with a
#'   warning. As with `"simple"` and `"OADP"`, imputed genotypes are used only
#'   if the `gen_tibble` is set to use them (see [gt_set_imputed()]).
#' * `"simple"`: the product of the standardised genotypes, with missing
#'   genotypes set to 0 (the mean), by the loadings. An individual with a
#'   proportion f of missing genotypes is shrunk towards 0 by a factor of
#'   about 1 - f (e.g. to a fifth of its distance with 80% missing
#'   genotypes), so it should only be used with complete or imputed data.
#' * `"plink"`: the simple projection divided by the proportion of loci
#'   genotyped in each individual, which is how PLINK 1.9 projects the
#'   individuals outside the clusters of `--pca-cluster-names` (or
#'   `--pca-clusters`). It removes most of the shrinkage of `"simple"`. PLINK
#'   standardises the genotypes with the allele frequencies of all the
#'   individuals in its run, reference and projected, whereas `predict()` uses
#'   those of the PCA (see [gt_pca()] to choose them): the scores are the same
#'   as PLINK's when the PCA uses the same frequencies, and PLINK's
#'   `.eigenvec` values are these scores divided by the singular values
#'   (`object$d`). Individuals with no genotyped locus get `NA`.
#' * `"OADP"`: the Online Augmentation, Decomposition, and Procrustes
#'   projection (Zhang et al. 2020), which corrects the shrinkage of
#'   projections with many more loci than individuals; missing genotypes are
#'   set to 0, as for `"simple"`, with the same shrinkage.
#' * `"none"`: the product of the genotypes by the loadings with
#'   [bigstatsr::big_prodMat()]; it needs complete or imputed genotypes (and
#'   uses imputed genotypes automatically).
#'
#' In a simulation of individuals with 20-80% missing genotypes projected on
#' complete reference individuals, `"least_squares"` had no shrinkage and the
#' smallest error, `"simple"` and `"OADP"` shrank individuals by their
#' proportion of missing genotypes, and imputing the genotypes (by their mode)
#' before a simple projection only halved the shrinkage. `"plink"` was close
#' to `"least_squares"` with the frequencies of the reference individuals, and
#' had a larger error with those of all the individuals, as PLINK does.
#'
#' @param object the [`gt_pca`] object
#' @param new_data a gen_tibble if scores are requested for a new dataset
#' @param project_method the projection method: `"least_squares"` (the
#'   default, as done by smartpca), `"simple"`, `"plink"` (as PLINK 1.9),
#'   `"OADP"` (Online Augmentation, Decomposition, and Procrustes projection)
#'   or `"none"`; see Details.
#' @param lsq_pcs a vector of the indices of the principal components to use
#'   for the least squares fitting (and to return). Only relevant if
#'   `project_method = "least_squares"`. It defaults to all the components of
#'   the PCA.
#' @param block_size number of loci read simultaneously (larger values will
#' speed up computation, but require more memory); by default, blocks of about
#' 128 MB.
#' @param n_cores number of cores
#' @param as_matrix logical, whether to return the result as a matrix (default)
#'   or a tibble.
#' @param ... no used
#' @returns a matrix of predictions (in line with predict using a prcomp
#'   object) or a tibble, with samples as rows and components as columns. The
#'   number of components depends on how many were estimated in the [`gt_pca`]
#'   object.
#' @references Zhang et al (2020). Fast and robust ancestry prediction using
#'   principal component analysis. Bioinformatics 36(11): 3439–3446.
#'
#'   Patterson N., Price A.L. and Reich D. (2006) Population structure and
#'   eigenanalysis. PLoS Genetics 2: e190.
#'
#'   Privé F., Luu K., Blum M.G.B., McGrath J.J. and Vilhjálmsson B.J. (2020)
#'   Efficient toolkit implementing best practices for principal component
#'   analysis of population genetic data. Bioinformatics 36: 4449-4457 (for
#'   the implementation of OADP in bigsnpr, which `"OADP"` uses).
#' @rdname predict_gt_pca
#' @importFrom foreach %do%
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' lobsters <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Remove monomorphic loci and impute
#' lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
#' lobsters <- gt_impute_simple(lobsters, method = "mode")
#'
#' # Subset into two datasets: one original and one to predict
#' original_lobsters <- lobsters[c(1:150), ]
#' new_lobsters <- lobsters[c(151:176), ]
#'
#' # Create PCA object
#' pca <- gt_pca(original_lobsters)
#'
#' # Predict
#' predict(pca, new_data = new_lobsters, project_method = "simple")
#'
#' # Predict as PLINK 1.9
#' predict(pca, new_data = new_lobsters, project_method = "plink")
#'
#' # Predict with OADP
#' predict(pca, new_data = new_lobsters, project_method = "OADP")
#'
#' # Predict with least squares (the default), on the first three components
#' predict(pca, new_data = new_lobsters, lsq_pcs = c(1, 2, 3))
#'
#' # Return a tibble
#' predict(pca, new_data = new_lobsters, as_matrix = FALSE)
#'
#' # Adjust block.size
#' predict(pca, new_data = new_lobsters, block_size = 10)
#'
# this is a modified version of bigstatsr::predict.big_SVD
predict.gt_pca <- function(
  object,
  new_data = NULL,
  project_method = c(
    "least_squares",
    "simple",
    "plink",
    "OADP",
    "none"
  ),
  lsq_pcs = NULL,
  block_size = NULL,
  n_cores = 1,
  as_matrix = TRUE,
  ...
) {
  if (n_cores > 1) {
    # Remove checking for two levels of parallelism
    .old_opt <- getOption("bigstatsr.check.parallel.blas", TRUE)
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = .old_opt), add = TRUE)
  }

  rlang::check_dots_empty()
  project_method <- match.arg(project_method)

  if (is.null(new_data)) {
    # U * D #nolint start
    UD <- sweep(object$u, 2, object$d, "*")
    dimnames(UD) <- list(rownames(object$u), paste0(".PC", seq_len(ncol(UD)))) # nolint end
    output <- output_type(object = UD, as_matrix, id = rownames(object$u))
    return(output)
  } else {
    if (!inherits(new_data, "gen_tbl")) {
      stop("new_data should be a gen_tibble")
    }
    # the FBMs of the projections are temporary (the result is a matrix or a
    # tibble): make them in a private directory, deleted on exit, so that
    # their names cannot clash with other .bk files (see
    # with_private_fbm_dir())
    local_private_fbm_dir()
    # get id of loci in new_data, and check that new_data have all the loci
    # used to build the pca (a single match() hashes the loci names once)
    loci_subset <- match(object$loci$name, show_loci(new_data)$name)
    if (anyNA(loci_subset)) {
      stop("loci used in object are not present in new_data")
    }
    if (
      !all(
        all(
          show_loci(new_data)$allele_ref[loci_subset] == object$loci$allele_ref
        ), # nolint
        all(
          show_loci(new_data)$allele_alt[loci_subset] == object$loci$allele_alt
        )
      )
    ) {
      stop(paste(
        "ref and alt alleles differ between new_data and the data",
        "used to create the pca object"
      ))
    }

    if (project_method == "none") {
      if (gt_has_imputed(new_data) && !gt_uses_imputed(new_data)) {
        gt_set_imputed(new_data, set = TRUE)
        on.exit(gt_set_imputed(new_data, set = FALSE), add = TRUE)
      }

      if (is.null(block_size)) {
        block_size <- bigstatsr::block_size(nrow(new_data))
      }
      # loci with a scale of 0 (monomorphic loci, which gt_pca() keeps, as
      # PLINK does) have loadings of 0: a scale of 1 avoids 0 / 0 (NaN)
      scale <- object$scale
      scale[scale == 0] <- 1
      # X * V #nolint
      XV <- bigstatsr::big_prodMat(
        .gt_get_fbm(new_data), # nolint
        object$v,
        ind.row = .gt_fbm_rows(new_data),
        ind.col = .gt_fbm_cols(new_data)[loci_subset],
        block.size = block_size,
        center = object$center,
        scale = scale
      )
      dimnames(XV) <- list(new_data$id, paste0(".PC", seq_len(ncol(XV)))) # nolint
      output <- output_type(object = XV, as_matrix, id = new_data$id)
      return(output)
    } else {
      # "simple", "plink", "OADP" and "least_squares": one pass over the loci of
      # new_data, in blocks (see pca_project_cpp())
      if (project_method == "least_squares") {
        # all the components by default
        if (is.null(lsq_pcs)) lsq_pcs <- seq_len(ncol(object$v))
        if (
          length(lsq_pcs) == 0 ||
            anyNA(lsq_pcs) ||
            any(lsq_pcs < 1) ||
            any(lsq_pcs > ncol(object$v)) ||
            !all(lsq_pcs == as.integer(lsq_pcs))
        ) {
          stop(
            "lsq_pcs should be a vector of valid component indices ",
            "(positive integers between 1 and ", ncol(object$v), "), ",
            "e.g., c(1, 2) or c(1, 2, 3)"
          )
        }
        # test that components are unique (i.e. no duplicates)
        if (any(duplicated(lsq_pcs))) {
          stop("lsq_pcs should not contain duplicate values")
        }
      }
      proj <- pca_project(object, new_data, loci_subset,
        lsq_pcs = if (project_method == "least_squares") {
          lsq_pcs
        } else {
          integer(0)
        },
        block_size = block_size, n_cores = n_cores
      )
      if (project_method == "simple") {
        XV <- proj$XV # nolint
      } else if (project_method == "plink") {
        # PLINK 1.9 divides the simple projection of each individual by its
        # number of genotyped loci (rather than the number of loci of the
        # PCA): divide by the proportion of loci genotyped (0 gives NA)
        obs_frac <- 1 - proj$n_missing / nrow(object$v)
        obs_frac[obs_frac == 0] <- NA
        XV <- proj$XV / obs_frac # nolint
      } else if (project_method == "OADP") {
        XV <- utils::getFromNamespace("OADP_proj", "bigsnpr")( # nolint
          proj$XV,
          proj$row_sums_sq,
          object$d,
          ncores = n_cores
        )
      } else {
        # For each individual, we solve the least squares problem
        #   crossprod(v_sub) b = crossprod(v_sub, genotypes_scaled)
        # where v_sub are the loadings (for lsq_pcs) of the loci genotyped in
        # that individual: crossprod(v_sub, genotypes_scaled) is the simple
        # projection (where missing genotypes contribute 0), and
        # crossprod(v_sub) is crossprod(V) minus the sum of v_j v_j' over the
        # loci j missing in that individual (proj$miss)
        n_lsq <- length(lsq_pcs)
        vtv <- crossprod(unclass(object$v)[, lsq_pcs, drop = FALSE])
        xv_lsq <- proj$XV[, lsq_pcs, drop = FALSE]
        # individuals genotyped at fewer loci than there are components have
        # a singular system (computed as a difference, it can look invertible)
        n_obs <- nrow(object$v) - proj$n_missing
        lsq_proj <- vapply(
          seq_len(nrow(new_data)),
          function(i) {
            if (n_obs[i] < n_lsq) {
              return(rep(NA_real_, n_lsq))
            }
            # a nearly singular system (e.g. loci genotyped only where the
            # loadings are 0) is an error, with a tolerance above rounding
            tryCatch(
              solve(
                vtv - matrix(proj$miss[i, ], n_lsq, n_lsq), xv_lsq[i, ],
                tol = 1e-10
              ),
              error = function(e) rep(NA_real_, n_lsq)
            )
          },
          numeric(n_lsq)
        )
        # one row per individual (vapply gives a vector for one component)
        lsq_proj <- matrix(lsq_proj, nrow = nrow(new_data), byrow = TRUE)
        n_na <- sum(is.na(lsq_proj[, 1]))
        if (n_na > 0) {
          warning(
            n_na, " individuals have too few genotyped loci to fit ",
            n_lsq, " components: their scores are NA"
          )
        }
        XV <- lsq_proj # nolint
      }
      pcs <- if (project_method == "least_squares") {
        lsq_pcs
      } else {
        seq_len(ncol(XV))
      }
      dimnames(XV) <- list(new_data$id, paste0(".PC", pcs)) # nolint
      output <- output_type(object = XV, as_matrix, id = new_data$id)
      return(output)
    }
  }
}

# The projection of new_data on the loadings of a gt_pca (for the loci
# loci_subset of new_data, in the order of the loci of the pca), with the
# blocked kernel pca_project_cpp(): returns XV, the row sums of squares (for
# OADP), the number of missing genotypes of each individual, and, for the
# components lsq_pcs, the crossproducts of the loadings of the loci missing in
# each individual (for least squares)
pca_project <- function(object, new_data, loci_subset, lsq_pcs, block_size,
                        n_cores) {
  v <- unclass(object$v)
  attributes(v) <- list(dim = dim(v))
  if (is.null(block_size)) block_size <- grm_block_size(nrow(new_data))
  pca_project_cpp(
    BM = .gt_get_fbm(new_data),
    rowInd = .gt_fbm_rows(new_data),
    colInd = .gt_fbm_cols(new_data)[loci_subset],
    center = object$center,
    scale = object$scale,
    V = v,
    lsq_pcs = as.integer(lsq_pcs),
    block_size = block_size,
    use_blas = grm_use_blas("auto", n_cores),
    ncores = n_cores
  )
}

# function to format the matrix into a tibble (if requested with
# as_matrix = FALSE)
output_type <- function(object, as_matrix, id) {
  if (!as_matrix) {
    object <- as_tibble(object)
    object <- object %>%
      dplyr::mutate(object, id = id) %>%
      dplyr::relocate(id)
    return(object)
  } else {
    return(object)
  }
}

###############################################################################
