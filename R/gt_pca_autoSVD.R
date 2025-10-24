#' PCA controlling for LD for `gen_tibble` objects
#'
#' This function performs Principal Component Analysis on a `gen_tibble`, using
#' a fast truncated SVD with initial pruning and then iterative removal of
#' long-range LD regions. This function is a wrapper for
#' [bigsnpr::snp_autoSVD()]
#'
#' Using gt_pca_autoSVD requires a reasonably large dataset, as the function
#' iteratively removes regions of long range LD. If you encounter: 'Error in
#' rollmean(): Parameter 'size' is too large.', `roll_size` exceeds the number
#' of variants on at least one of your chromosomes. Try reducing 'roll_size' to
#' avoid this error.
#'
#' Note: rather than accessing these elements directly, it is better to use
#' `tidy` and `augment`. See [`gt_pca_tidiers`].
#'
#' @param x a `gen_tbl` object
#' @param k Number of singular vectors/values to compute. Default is `10`.
#'   **This algorithm should be used to compute a few singular vectors/values.**
#' @param fun_scaling Usually this  can be left unset, as it defaults to
#'   [bigsnpr::snp_scaleBinom()], which is the appropriate function for
#'   biallelic SNPs. Alternatively it is possible to use  custom function (see
#'   [bigsnpr::snp_autoSVD()] for details.
#' @param total_var a boolean indicating whether to compute the total variance
#'   of the matrix. Default is `TRUE`. Using `FALSE` will speed up computation,
#'   but the total variance will not be stored in the output (and thus it will
#'   not be possible to assign a proportion of variance explained to the
#'   components).
#' @param roll_size Radius of rolling windows to smooth log-p-values. Default is
#'   `50`.
#' @param int_min_size Minimum number of consecutive outlier SNPs in order to be
#'   reported as long-range LD region. Default is `20`.
#' @param verbose Output some information on the iterations? Default is `TRUE`.
#' @param thr_r2 Threshold over the squared correlation between two SNPs.
#'   Default is `0.2`. Use `NA` if you want to skip the clumping step. size
#' @param use_positions a boolean on whether the position is used to define
#'   `size`, or whether the size should be in number of SNPs. Default is TRUE
#' @param size For one SNP, window size around this SNP to compute correlations.
#'   Default is 100 / thr_r2 for clumping (0.2 -> 500; 0.1 -> 1000; 0.5 -> 200).
#'   If not providing infos.pos (NULL, the default), this is a window in number
#'   of SNPs, otherwise it is a window in kb (genetic distance). I recommend
#'   that you provide the positions if available.
#' @param alpha_tukey Default is `0.05`. The type-I error rate in outlier
#'   detection (that is further corrected for multiple testing).
#' @param min_mac Minimum minor allele count (MAC) for variants to be included.
#'   Default is `10`.
#' @param max_iter Maximum number of iterations of outlier detection. Default is
#'   `5`.
#' @param n_cores Number of cores used. Default doesn't use parallelism. You may
#'   use [bigstatsr::nb_cores()].
#' @returns a `gt_pca` object, which is a subclass of `bigSVD`; this is an S3
#'   list with elements: A named list (an S3 class "big_SVD") of
#' - `d`, the eigenvalues (singular values, i.e. as variances),
#' - `u`, the scores for each sample on each component
#'   (the left singular vectors)
#' - `v`, the loadings (the right singular vectors)
#' - `center`, the centering vector,
#' - `scale`, the scaling vector,
#' - `method`, a string defining the method (in this case 'autoSVD'),
#' - `call`, the call that generated the object.
#' - `loci`, the loci used after long range LD removal.
#' @seealso [bigsnpr::snp_autoSVD()]  which this function wraps.
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
#' show_loci(lobsters)$chromosome <- "1"
#'
#' # Create PCA object, including total variance
#' gt_pca_autoSVD(lobsters,
#'   k = 10,
#'   roll_size = 20,
#'   total_var = TRUE
#' )
#' # Change number of components and exclude total variance
#' gt_pca_autoSVD(lobsters,
#'   k = 5,
#'   roll_size = 20,
#'   total_var = FALSE
#' )
#'
# nolint start
gt_pca_autoSVD <- function(
    # nolint end
    x,
    k = 10,
    fun_scaling = bigsnpr::snp_scaleBinom(),
    thr_r2 = 0.2,
    use_positions = TRUE,
    size = 100 / thr_r2,
    roll_size = 50,
    int_min_size = 20,
    alpha_tukey = 0.05,
    min_mac = 10,
    max_iter = 5,
    n_cores = 1,
    verbose = TRUE,
    total_var = TRUE) {
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

  X <- attr(x$genotypes, "fbm") # convenient pointer #nolint
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  # we need to create a chromosome vector that is as long as
  # the complete bigsnp object
  infos_chr <- rep(1, ncol(X))

  infos_chr[.gt_fbm_cols(x)] <- cast_chromosome_to_int(show_loci(x)$chromosome)
  # chromosomes have to be positive numbers
  if (min(infos_chr) < 1) {
    infos_chr <- infos_chr + abs(min(infos_chr) + 1)
  }
  infos_pos <- NULL
  if (use_positions) {
    is_loci_table_ordered(x, error_on_false = TRUE)
    infos_pos <- rep(0, ncol(X))
    infos_pos[.gt_fbm_cols(x)] <- show_loci(x)$position
  }
  # TODO
  # Do we want to use the code from loci_clump to create chromosomes and
  # positions (it is a bit neater)

  tryCatch(
    expr = {
      this_svd <- bigsnpr::snp_autoSVD(
        X, # nolint
        infos.chr = infos_chr,
        infos.pos = infos_pos,
        ind.row = .gt_fbm_rows(x),
        ind.col = .gt_fbm_cols(x),
        fun.scaling = fun_scaling,
        thr.r2 = thr_r2,
        size = size,
        k = k,
        roll.size = roll_size,
        int.min.size = int_min_size,
        alpha.tukey = alpha_tukey,
        min.mac = min_mac,
        max.iter = max_iter,
        ncores = n_cores,
        verbose = verbose
      )
    },
    error = function(e) {
      if (grepl("Parameter 'size' is too large.",
        e$message,
        fixed = TRUE
      )) {
        stop(
          "'Error in rollmean(): Parameter 'size' is too large.'
          roll_size exceeds the number of variants on at least one of your
          chromosomes. Try reducing 'roll_size' to avoid this error. "
        )
      } else {
        stop(e)
      }
    }
  )

  # add names to the scores (to match them to data later)
  rownames(this_svd$u) <- x$id
  this_svd$method <- "autoSVD"
  this_svd$call <- match.call()
  # subset the loci table to have only the snps of interest
  this_svd$loci <-
    show_loci(x)[.gt_fbm_cols(x) %in% attr(this_svd, "subset"), ]
  class(this_svd) <- c("gt_pca", class(this_svd))
  if (total_var) {
    loci <- this_svd$loci
    loci_after_ld <- which(show_loci(x)$name %in% loci$name)
    x_autoSVD_subset <- x %>% select_loci(all_of(loci_after_ld)) # nolint
    x_ind_col <- show_loci(x_autoSVD_subset)$big_index
    x_ind_row <- vctrs::vec_data(x_autoSVD_subset$genotypes)
    this_svd$square_frobenius <- square_frobenius(
      X, # nolint
      x_ind_row,
      x_ind_col,
      center = this_svd$center,
      scale = this_svd$scale
    )
  }

  this_svd
}
