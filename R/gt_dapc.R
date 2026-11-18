#' Discriminant Analysis of Principal Components for gen_tibble
#'
#' This function implements the Discriminant Analysis of Principal Components
#' (DAPC, Jombart et al. 2010). This method describes the diversity between
#' pre-defined groups. When groups are unknown, use [gt_cluster_pca()] to infer
#' genetic clusters. See 'details' section for a succinct description of the
#' method, and the vignette in the package `adegenet` ("adegenet-dapc") for a
#' tutorial.
#'
#' The Discriminant Analysis of Principal Components (DAPC) is designed to
#' investigate the genetic structure of biological populations. This
#' multivariate method consists in a two-steps procedure. First, genetic data
#' are transformed (centred, possibly scaled) and submitted to a Principal
#' Component Analysis (PCA). Second, principal components of PCA are submitted
#' to a Linear Discriminant Analysis (LDA). A trivial matrix operation allows to
#' express discriminant functions as linear combination of alleles, therefore
#' allowing one to compute allele contributions. More details about the
#' computation of DAPC are to be found in the indicated reference.
#'
#' Results can be visualised with [`autoplot.gt_dapc()`], see the help of that
#' method for the available plots. There are also [gt_dapc_tidiers] for
#' manipulating the results. For the moment, this function returns objects of
#' class [`adegenet::dapc`] which are
#' compatible with methods from `adegenet`; graphical methods for DAPC are
#' documented in [adegenet::scatter.dapc] (see ?scatter.dapc). This is likely
#' to change in the future, so make sure you do not rely on the objects
#' remaining compatible.
#'
#' This function aligns with the guidelines proposed by Thia (2023) for the
#' standardized application of DAPC to genotype data. Our default settings are
#' designed to follow these recommendations, so that the number of principal
#' components (`n_pca`) defaults to the smaller of *k*-1 and the number of
#' available principal components (where *k* is the number of populations or
#' clusters), and the number of discriminant functions (`n_da`) is set to the
#' minimum of *k*-1 and `n_pca`. The user can override these defaults by
#' specifying the `n_pca` and `n_da` arguments, but caution is advised when
#' adjusting `n_pca` to avoid potential overfitting. We recommend users consult
#' these guidelines and consider their individual dataset to ensure best
#' practices.
#'
#' With missing genotypes, the PCA does not need imputed genotypes: the
#' default method of [gt_pca()] compares individuals over the loci genotyped
#' in both (as PLINK does), and the discriminant analysis uses its scores as
#' it would those of any PCA. The loadings of the loci (`loadings_by_locus`)
#' are then estimated by regressing the standardised genotypes of each locus
#' on the scores of the `n_pca` components (with an intercept), over the
#' individuals genotyped at that locus, and multiplying these coefficients by
#' the loadings of the discriminant functions. This needs the genotypes: give
#' the `gen_tibble` of the PCA as `data`. With complete genotypes, the
#' regression coefficients are the loadings of the PCA, so the results are
#' the same as with the usual matrix product. With missing genotypes, the
#' loadings of the PCA (which count missing genotypes as 0) would shrink loci
#' with many missing genotypes, and hide them in the loadings plots. Loci
#' genotyped in fewer than `n_pca + 1` individuals get `NA` loadings. The
#' posterior probabilities of group membership do not account for the
#' noisier scores of individuals with many missing genotypes.
#'
#' Note that there is no current method to predict scores for
#' individuals not included in the original analysis. This is because we
#' currently do not have a mechanism to store the pca information in the
#' object, and that is needed for prediction.
#'
#' @references Jombart T, Devillard S and Balloux F (2010) Discriminant analysis
#'   of principal components: a new method for the analysis of genetically
#'   structured populations. BMC Genetics 11:94. doi:10.1186/1471-2156-11-94
#'   Thia, J. A. (2023). Guidelines for standardizing the application of
#'   discriminant analysis of principal components to genotype data. Molecular
#'   Ecology Resources, 23, 523–538. https://doi.org/10.1111/1755-0998.13706
#'
#'
#' @param x an object of class `gt_pca`, or its subclass `gt_cluster_pca`
#' @param pop either a factor indicating the group membership of individuals; or
#'   an integer defining the desired *k* if x is a `gt_cluster_pca`; or NULL, if
#'   'x' is a `gt_cluster_pca` and contain an element 'best_k', usually
#'   generated with [gt_cluster_pca_best_k()], which will be used to select the
#'   clustering level.
#' @param n_pca number of principal components to be used in the Discriminant
#'   Analysis. If NULL, k-1 will be used.
#' @param n_da an integer indicating the number of axes retained in the
#'   Discriminant Analysis step.
#' @param loadings_by_locus a logical indicating whether the loadings and
#'   contribution of each locus should be stored (TRUE, default) or not (FALSE).
#'   Such output can be useful, but can also create large matrices when there
#'   are a lot of loci and many dimensions.
#' @param data the `gen_tibble` the PCA `x` was computed on (it can have more
#'   loci than the PCA). It is needed only for the loadings of the loci of a
#'   PCA with missing genotypes (see Details), and ignored otherwise.
#' @param n_cores the number of cores used to estimate the loadings of the
#'   loci of a PCA with missing genotypes. Default is
#'   [bigstatsr::nb_cores()].
#' @returns an object of class [adegenet::dapc]
#' @export
#' @seealso [gt_cluster_pca()] [gt_cluster_pca_best_k()] [adegenet::dapc()]
#' @examples
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' lobsters <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Remove monomorphic loci
#' lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
#'
#' # Create PCA object (the lobsters have some missing genotypes, which do not
#' # need to be imputed)
#' pca <- gt_pca(lobsters)
#'
#' # Run DAPC on the `gt_pca` object, providing `pop` as factor; with missing
#' # genotypes, the loadings of the loci need the genotypes (`data`)
#' populations <- as.factor(lobsters$population)
#' gt_dapc(pca, n_pca = 6, n_da = 2, pop = populations, data = lobsters)
#'
#' # Run clustering on the first 10 PCs
#' cluster_pca <- gt_cluster_pca(
#'   x = pca,
#'   n_pca = 10,
#'   k_clusters = c(1, 5),
#'   method = "kmeans",
#'   n_iter = 1e5,
#'   n_start = 10,
#'   quiet = FALSE
#' )
#'
#' # Find best k
#' cluster_pca <- gt_cluster_pca_best_k(cluster_pca,
#'   stat = "BIC",
#'   criterion = "min"
#' )
#'
#' # Run DAPC on the `gt_cluster_pca` object
#' gt_dapc(cluster_pca, n_pca = 10, n_da = 2, data = lobsters)
#'
# @param pca_info a logical indicating whether information about the prior PCA
#' #  should be stored (TRUE, default) or not (FALSE). This information is
#' #  required to predict group membership of new individuals using predict, but
#' #  makes the object slightly bigger.
# AM: Thoughts about data structures. The original DAPC blended pca info
# within the object. For a cleaner job at predicting, it would be best to
# simply store the pca object as an element within the object. This would break
# compatibility, but only with the predict function, which in any case will not
# work.
# Once we have tidiers and autoplot, it might be best to reorganise the
# object to fully fit our purposes, and give up trying to be backcompatible
# with adegenet
gt_dapc <- function(
  x,
  pop = NULL,
  n_pca = NULL,
  n_da = NULL,
  loadings_by_locus = TRUE,
  data = NULL,
  n_cores = bigstatsr::nb_cores()
) {
  if (!inherits(x, "gt_pca")) {
    stop("'x' should be a 'gt_pca' object")
  }
  if (is.null(x$center)) {
    stop("'x' was run without centering; centering is necessary for 'gt_dapc'")
  }
  # with a PCA of missing genotypes, the loadings of the loci are estimated
  # from the genotypes (see loci_lsq_loadings()), so they are needed
  pca_missing <- isTRUE(attr(x, "missing"))
  if (loadings_by_locus && pca_missing && is.null(data)) {
    stop(
      "'x' is a PCA of genotypes with missing values: the loadings of the ",
      "loci need the genotypes, so give the gen_tibble of the PCA as ",
      "'data', or set 'loadings_by_locus = FALSE'"
    )
  }

  pop.fac <- NULL
  if (is.null(pop)) {
    # if no pop was given, use best_k
    if (any(!inherits(x, "gt_cluster_pca"), is.null(x$best_k))) {
      stop("if 'pop' is not set, 'x' should be a 'gt_cluster_pca'")
    }
    # groups are named by their k value
    pop.fac <- x$clusters$groups[[as.character(x$best_k)]]
  } else if (is.factor(pop)) {
    # if a factor with all assignments was given
    pop.fac <- pop
  } else if (is.numeric(pop) && inherits(x, "gt_cluster_pca")) {
    # if pop is the k value
    pop.fac <- x$clusters$groups[[as.character(pop)]]
  }

  if (is.null(pop.fac)) {
    stop(paste(
      "'pop' should be a factor, or a value of k in 'x$clusters$k'",
      "if 'x' is a 'gt_cluster_pca'"
    ))
  }
  pop.fac <- as.factor(pop.fac)
  # FIXME: is this actually nlevels(pop)?
  n_pop <- nlevels(pop.fac)

  if (is.null(n_pca)) {
    # if we generated clusters, use the same pca
    if (inherits(x, "gt_cluster_pca")) {
      n_pca <- x$clusters$n_pca
    } else {
      # use all principal components
      n_pca <- length(x$d)
    }
    # by default, use number of clusters minus 1
    if (n_pca > n_pop - 1) {
      n_pca <- n_pop - 1
    }
  } else {
    # if n_pca was given, check that it is not too large
    if (n_pca > ncol(x$u)) {
      n_pca <- ncol(x$u)
    }
  }

  XU <- sweep(x$u, 2, x$d, "*")[, 1:n_pca, drop = FALSE] # principal components
  # note that this is the proportion of variance out of the variance
  # we started with (i.e. what we retained with the PCAs)
  # sum of retained eigenvalues (the squared singular values)
  XU.lambda <- sum(x$d[1:n_pca]^2) / sum(x$d^2)
  colnames(XU) <- paste("PCA-pc", seq_len(ncol(XU)), sep = ".")

  ## PERFORM DA ##
  # tol=1e-30 is a kludge, but a safe (?) one to avoid fancy
  # rescaling by lda.default
  ldaX <- MASS::lda(XU, pop.fac, tol = 1e-30)
  lda.dim <- sum(ldaX$svd^2 > 1e-10)
  ldaX$svd <- ldaX$svd[1:lda.dim]
  ldaX$scaling <- ldaX$scaling[, 1:lda.dim, drop = FALSE]

  # can't be more than K-1 disc. func., or more than n.pca
  if (is.null(n_da)) {
    n_da <- min(length(levels(pop.fac)) - 1, n_pca, sum(ldaX$svd > 1e-10))
  } else {
    n_da <- min(n_da, length(levels(pop.fac)) - 1, n_pca, sum(ldaX$svd > 1e-10))
  }
  n_da <- round(n_da)
  predX <- stats::predict(ldaX, dimen = n_da)

  ## BUILD RESULT
  res <- list()
  res$n.pca <- n_pca
  res$n.da <- n_da
  res$tab <- XU
  res$grp <- pop.fac
  # note that this (res$var) is the variance out of the variance that was
  # captured by the retained PCs
  res$var <- XU.lambda
  res$eig <- ldaX$svd^2
  res$loadings <- ldaX$scaling[, 1:n_da, drop = FALSE]
  res$means <- ldaX$means
  res$ind.coord <- predX$x
  res$grp.coord <- apply(res$ind.coord, 2, tapply, pop.fac, mean)
  res$prior <- ldaX$prior
  res$posterior <- predX$posterior
  res$assign <- predX$class
  res$call <- match.call()

  # nolint start
  # ## optional: store loadings of variables
  # @BUG we need to sort out the slots as these are not correct
  # @TODO our objects are missing several of these slots
  # if (pca_info) {
  #   warning(paste(
  #     "conversion of objects slots is inconmplete, don't use",
  #     "this option yet!"
  #   ))
  #   res$pca.loadings <- as.matrix(V)
  #   # res$pca.cent <- x$cent
  #   # if(!is.null(x$norm)) {
  #   #   res$pca.norm <- x$norm
  #   # } else {
  #   #   res$pca.norm <- rep(1, length(x$cent))
  #   # }
  #   res$pca.eig <- x$d^2 # TODO check, this should get back the eigen from glPCA
  #   # note that the default allele.as.unit is FALSE for glPCA
  # }
  # nolint end

  ## optional: get loadings of variables
  if (loadings_by_locus) {
    V <- x$v[, 1:n_pca, drop = FALSE] # principal axes
    if (pca_missing) {
      # with missing genotypes, the loadings of the PCA (missing genotypes as
      # 0) are biased: regress each locus on the scores instead, over the
      # individuals genotyped at that locus (with complete genotypes, this
      # gives the loadings of the PCA)
      V[] <- loci_lsq_loadings(x, data, n_pca = n_pca, n_cores = n_cores)
    }
    colnames(V) <- paste("PCA-pa", seq_len(ncol(V)), sep = ".")
    var.load <- as.matrix(V) %*% as.matrix(ldaX$scaling[, 1:n_da, drop = FALSE])
    f1 <- function(x) {
      # (loci whose loadings could not be estimated are NA, and are left out
      # of the sum)
      temp <- sum(x * x, na.rm = TRUE)
      if (temp < 1e-12) {
        return(rep(0, length(x)))
      }
      return(x * x / temp)
    }
    # matrix multiplication and apply() drop the integer loci names of a
    # matrix_int_names, so we add them back
    res$var.contr <- matrix_int_names(apply(var.load, 2, f1),
      row_names = row_names(V)
    )
    res$var.load <- matrix_int_names(var.load, row_names = row_names(V))
  }

  class(res) <- c("gt_dapc", "dapc")
  return(res)
}

# Loadings of the loci for gt_dapc() with a PCA of missing genotypes: the
# coefficients of the least squares regression of the standardised genotypes
# of each locus (with the centre and scale of the PCA) on an intercept and
# the scores (u d) of the first n_pca components, over the individuals
# genotyped at that locus (see grm_loci_lsq_cpp()). With complete genotypes,
# these are the loadings v of the PCA; with missing genotypes, the loadings
# of the PCA (missing genotypes as 0) shrink each locus by the fraction of
# individuals genotyped there, and are distorted when missingness differs
# between groups.
# @param pca a gt_pca object
# @param data the gen_tibble the PCA was computed on (it can have more loci)
# @param n_pca the number of components
# @param n_cores the number of cores
# @returns a matrix with one row per locus of the PCA and n_pca columns, NA
#   for loci genotyped in too few individuals (with a warning)
# @keywords internal
# @noRd
loci_lsq_loadings <- function(pca, data, n_pca, n_cores) {
  if (!inherits(data, "gen_tbl")) {
    stop("'data' should be a gen_tibble")
  }
  # the individuals must be those of the PCA, in the same order
  if (!identical(as.character(data$id), rownames(pca$u))) {
    stop("the individuals in 'data' do not correspond to those in the PCA")
  }
  # the loci of the PCA in data, with the same alleles
  loci_subset <- pca_match_loci(pca, data,
    what = "'data'", pca_name = "the PCA"
  )
  # the PCA used the genotypes with their missing values, not imputed ones
  if (gt_has_imputed(data) && gt_uses_imputed(data)) {
    gt_set_imputed(data, set = FALSE)
    on.exit(gt_set_imputed(data, set = TRUE), add = TRUE)
  }
  n_ind <- nrow(data)
  pcs <- seq_len(n_pca)
  # the regression is on u sqrt(n), whose columns have the same norm as the
  # intercept, so that the test of nearly singular systems (in
  # grm_loci_lsq_cpp()) is about the individuals genotyped at a locus, not
  # about the scale of the scores
  u <- unclass(pca$u)[, pcs, drop = FALSE]
  coefs <- grm_loci_lsq_cpp(
    BM = .gt_get_fbm(data),
    rowInd = .gt_fbm_rows(data),
    colInd = .gt_fbm_cols(data)[loci_subset],
    center = pca$center,
    scale = pca$scale,
    S = u * sqrt(n_ind),
    block_size = grm_block_size(n_ind),
    ncores = n_cores,
    use_blas = grm_use_blas("auto", n_cores)
  )
  # rescale the coefficients to the scores u d:
  # z = a + (u sqrt(n)) c = a + (u d) (c sqrt(n) / d)
  coefs <- sweep(coefs, 2, sqrt(n_ind) / pca$d[pcs], "*")
  n_na <- sum(is.na(coefs[, 1]))
  if (n_na > 0) {
    warning(
      n_na, " loci are genotyped in too few individuals to estimate their ",
      "loadings (at least n_pca + 1 are needed, with varied scores): ",
      "their loadings are NA"
    )
  }
  coefs
}
