#' Predict scores of a PCA
#'
#' Predict the PCA scores for a [`gt_pca`], either for the original data or
#' projecting new data.
#'
#' @param object the [`gt_pca`] object
#' @param new_data a gen_tibble if scores are requested for a new dataset
#' @param project_method a string taking the value of either "simple", "OADP"
#'   (Online Augmentation, Decomposition, and Procrustes (OADP) projection), or
#'   "least_squares" (as done by SMARTPCA)
#' @param lsq_pcs a vector of the indices of the principal
#' components to use for the least square fitting. Only relevant
#' if `project_method = 'least_squares'`. It defaults to the first two
#' components.
#' @param block_size number of loci read simultaneously (larger values will
#' speed up computation, but require more memory)
#' @param n_cores number of cores
#' @param as_matrix logical, whether to return the result as a matrix (default)
#'   or a tibble.
#' @param ... no used
#' @returns a matrix of predictions (in line with predict using a prcomp
#'   object) or a tibble, with samples as rows and components as columns. The
#'   number of components depends on how many were estimated in the [`gt_pca`]
#'   object.
#' @references Zhang et al (2020). Fast and robust ancestry prediction using
#'   principal component analysis  36(11): 3439–3446.
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
#' # Predict with OADP
#' predict(pca, new_data = new_lobsters, project_method = "OADP")
#'
#' # Predict with least squares
#' predict(pca,
#'   new_data = new_lobsters,
#'   project_method = "least_squares", lsq_pcs = c(1, 2, 3)
#' )
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
    "none",
    "simple",
    "OADP",
    "least_squares"
  ),
  lsq_pcs = c(1, 2),
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
      # X * V #nolint
      XV <- bigstatsr::big_prodMat(
        .gt_get_fbm(new_data), # nolint
        object$v,
        ind.row = .gt_fbm_rows(new_data),
        ind.col = .gt_fbm_cols(new_data)[loci_subset],
        block.size = block_size,
        center = object$center,
        scale = object$scale
      )
      dimnames(XV) <- list(new_data$id, paste0(".PC", seq_len(ncol(XV)))) # nolint
      output <- output_type(object = XV, as_matrix, id = new_data$id)
      return(output)
    } else {
      # "simple", "OADP" and "least_squares": one pass over the loci of
      # new_data, in blocks (see pca_project_cpp())
      if (project_method == "least_squares") {
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
        lsq_proj <- vapply(
          seq_len(nrow(new_data)),
          function(i) {
            solve(vtv - matrix(proj$miss[i, ], n_lsq, n_lsq), xv_lsq[i, ])
          },
          numeric(n_lsq)
        )
        # one row per individual
        XV <- matrix(lsq_proj, nrow = nrow(new_data), byrow = TRUE) # nolint
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
