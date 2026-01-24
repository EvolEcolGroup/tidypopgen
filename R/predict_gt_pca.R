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
#' pca <- gt_pca_partialSVD(original_lobsters)
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
    # check the new_data have the same loci as the dataset used to build the pca
    if (!all(object$loci$name %in% show_loci(new_data)$name)) {
      stop("loci used in object are not present in new_data")
    }
    # get id of loci in new_data
    loci_subset <- match(object$loci$name, show_loci(new_data)$name)
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
    } else if (any(c("simple", "OADP") %in% project_method)) {
      X_norm <- bigstatsr::FBM(nrow(new_data), 1, init = 0) # nolint
      XV <- bigstatsr::FBM(nrow(new_data), ncol(object$u), init = 0) # nolint

      bigstatsr::big_parallelize(
        .gt_get_fbm(new_data),
        p.FUN = fbm256_part_prod,
        ind = seq_along(loci_subset),
        ncores = n_cores,
        ind.row = .gt_fbm_rows(new_data),
        # info_snp$`_NUM_ID_`[keep], #nolint
        ind.col = .gt_fbm_cols(new_data)[loci_subset],
        center = object$center,
        scale = object$scale,
        V = object$v,
        XV = XV,
        X_norm = X_norm
      )
      if (project_method == "simple") {
        XV <- XV[, , drop = FALSE] # nolint
        dimnames(XV) <- list(new_data$id, paste0(".PC", seq_len(ncol(XV)))) # nolint
        output <- output_type(object = XV, as_matrix, id = new_data$id)
        return(output)
      } else {
        XV <- utils::getFromNamespace("OADP_proj", "bigsnpr")(
          XV,
          X_norm,
          object$d,
          ncores = n_cores
        ) # nolint
        dimnames(XV) <- list(new_data$id, paste0(".PC", seq_len(ncol(XV)))) # nolint
        output <- output_type(object = XV, as_matrix, id = new_data$id)
        return(output)
      }
    } else if (project_method == "least_squares") {
      if (
          length(lsq_pcs) == 0 ||
            anyNA(lsq_pcs) ||
            any(lsq_pcs < 1) ||
            any(lsq_pcs > ncol(object$v)) ||
            !all(lsq_pcs == as.integer(lsq_pcs))) {
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
      X <- .gt_get_fbm(new_data) # pointer for FBM #nolint
      proj_i <- NULL # hack to define iterator
      lsq_proj <- foreach::foreach(
        proj_i = seq_len(nrow(new_data)),
        .final = t,
        .combine = cbind
      ) %do% {
        # scaled genotypes for this individual
        genotypes <- X[
          .gt_fbm_rows(new_data)[proj_i],
          .gt_fbm_cols(new_data)[loci_subset]
        ]
        genotypes_scaled <- (genotypes - object$center) / object$scale
        na_ids <- which(!is.na(genotypes_scaled))
        genotypes_scaled <- genotypes_scaled[na_ids]
        v_sub <- object$v[na_ids, lsq_pcs]
        solve(crossprod(v_sub), crossprod(v_sub, genotypes_scaled))
      }
      dimnames(lsq_proj) <-
        list(new_data$id, paste0(".PC", lsq_pcs))
      output <- output_type(object = lsq_proj, as_matrix, id = new_data$id)
      return(output)
    }
  }
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
# a port of bigsnpr::part_prod to work on standard fb256 matrices

fbm256_part_prod <- function(
  X,
  ind,
  ind.row,
  ind.col,
  center, # nolint
  scale,
  V,
  XV,
  X_norm
) {
  # nolint
  res <- fbm256_prod_and_rowSumsSq(
    BM = X,
    ind_row = ind.row,
    ind_col = ind.col[ind],
    center = center[ind],
    scale = scale[ind],
    V = V[ind, , drop = FALSE]
  )

  bigstatsr::big_increment(XV, res[[1]], use_lock = TRUE)
  bigstatsr::big_increment(X_norm, res[[2]], use_lock = TRUE)
}
