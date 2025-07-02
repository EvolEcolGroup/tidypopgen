#' Tidy a `gt_pca` object
#'
#' This summarizes information about the components of a `gt_pca` from the
#' `tidypopgen` package. The parameter `matrix` determines which element is
#' returned. Column names of the tidied output match those returned by
#' [broom::tidy.prcomp], the tidier for the standard PCA objects returned
#' by [stats::prcomp].
#'
#' @param x A `gt_pca` object returned by one of the `gt_pca_*` functions.
#' @param matrix Character specifying which component of the PCA should be
#'   tidied.
#'
#'   - `"samples"`, `"scores"`, or `"x"`: returns information about
#'     the map from the original space into principle components space
#'     (this is equivalent to product of *u* and *d*).
#'
#'   - `"v"`, `"rotation"`, `"loadings"` or `"variables"`: returns information
#'     about the map from principle components space back into the original
#'     space.
#'
#'   - `"d"`, `"eigenvalues"` or `"pcs"`: returns information about the
#'     eigenvalues.
#'
#' @param ... Not used. Needed to match generic signature only.
#'
#' @return A [tibble::tibble] with columns depending on the component of
#'   PCA being tidied.
#'
#'   If `"scores"` each row in the
#'   tidied output corresponds to the original data in PCA space. The columns
#'   are:
#'
#'   \item{`row`}{ID of the original observation (i.e. rowname from original
#'     data).}
#'   \item{`PC`}{Integer indicating a principal component.}
#'   \item{`value`}{The score of the observation for that particular principal
#'     component. That is, the location of the observation in PCA space.}
#'
#'   If `matrix` is  `"loadings"`, each
#'   row in the tidied output corresponds to information about the principle
#'   components in the original space. The columns are:
#'
#'   \item{`row`}{The variable labels (colnames) of the data set on
#'   which PCA was performed.}
#'   \item{`PC`}{An integer vector indicating the principal component.}
#'   \item{`value`}{The value of the eigenvector (axis score) on the
#'   indicated principal component.}
#'
#'   If `"eigenvalues"`, the columns are:
#'
#'   \item{`PC`}{An integer vector indicating the principal component.}
#'   \item{`std.dev`}{Standard deviation (i.e. sqrt(eig/(n-1))) explained by
#'   this PC (for compatibility with `prcomp`.}
#'   \item{`cumulative`}{Cumulative variation explained by
#'     principal components up to this component (note that this is NOT
#'     phrased as a percentage of total variance, since many methods
#'     only estimate a truncated SVD.}
#'
#' @aliases gt_pca_tidiers
#' @rdname tidy_gt_pca
#' @export
#' @seealso [gt_pca_autoSVD()] [augment_gt_pca]
#' @examples
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
#' # Create PCA object
#' pca <- gt_pca_partialSVD(lobsters)
#'
#' # Tidy the PCA object
#' tidy(pca)
#'
#' # Tidy the PCA object for eigenvalues
#' tidy(pca, matrix = "eigenvalues")
#'
#' # Tidy the PCA object for loadings
#' tidy(pca, matrix = "loadings")
#'
#' # Tidy the PCA object for scores
#' tidy(pca, matrix = "scores")
tidy.gt_pca <- function(x, matrix = "eigenvalues", ...) {
  if (length(matrix) > 1) {
    stop("Must select a single matrix to tidy.", call. = FALSE)
  }

  matrix_args <- c(
    "rotation",
    "x",
    "variables",
    "samples",
    "v",
    "pcs",
    "d",
    "scores",
    "loadings",
    "eigenvalues"
  )
  matrix <- rlang::arg_match(matrix, matrix_args)

  if (matrix %in% c("pcs", "d", "eigenvalues")) {
    ret <- tibble(
      PC = seq_len(length(x$d)),
      "std.dev" = x$d / sqrt(nrow(x$u) - 1)
    )
    if ("square_frobenius" %in% names(x)) {
      rssq <- x$square_frobenius
      var_exp <- x$d^2 / rssq
      percentage <- var_exp * 100
      cum_percentage <- cumsum(var_exp) * 100
      ret <- ret %>%
        mutate(percent = percentage) %>%
        mutate(cumulative = cum_percentage)
    }
  } else if (matrix %in% c("rotation", "variables", "v", "loadings")) {
    colnames(x$v) <- seq_len(ncol(x$v))
    ret <- x$v %>%
      tibble::as_tibble(rownames = "column") %>%
      tidyr::pivot_longer(
        cols = -"column",
        names_to = "PC",
        values_to = "value"
      )
    ret <- mutate(ret, PC = as.numeric(.data$PC))
    if (is.null(rownames(x$v))) ret$column <- as.integer(ret$column)
  } else if (matrix %in% c("x", "samples", "scores")) {
    ret <- sweep(x$u, 2, x$d, "*")
    colnames(ret) <- paste0("PC", seq_len(ncol(ret)))
    ret <- ret %>%
      tibble::as_tibble(rownames = "row") %>%
      tidyr::pivot_longer(
        cols = -"row",
        names_to = "PC",
        values_to = "value"
      )
    if (is.null(rownames(x$u))) ret$row <- as.integer(ret$row)
  }

  # change the PC to a numeric
  ret <- mutate(ret, PC = as.integer(str_replace_base(.data$PC, "PC", "")))
  as_tibble(ret)
}


#' Augment data with information from a gt_pca object
#'
#' Augment for `gt_pca` accepts a model object and a dataset and adds scores to
#' each observation in the dataset. Scores for each component are stored in a
#' separate column, which is given name with the pattern ".fittedPC1",
#' ".fittedPC2", etc. For consistency with [broom::augment.prcomp], a column
#' ".rownames" is also returned; it is a copy of 'id', but it ensures that any
#' scripts written for data augmented with [broom::augment.prcomp] will work out
#' of the box (this is especially helpful when adapting plotting scripts).
#' @param x  A `gt_pca` object returned by one of the `gt_pca_*` functions.
#' @param data the `gen_tibble` used to run the PCA.
#' @param k the number of components to add
#' @param ... Not used. Needed to match generic signature only.
#' @return A  [gen_tibble] containing the original data along with additional
#'   columns containing each observation's projection into PCA space.
#' @export
#' @name augment_gt_pca
#' @seealso [gt_pca_autoSVD()] [gt_pca_tidiers]
#' @examples
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
#' # Create PCA object
#' pca <- gt_pca_partialSVD(lobsters)
#'
#' # Augment the gen_tibble with PCA scores
#' augment(pca, data = lobsters)
#'
#' # Adjust the number of components to add
#' augment(pca, data = lobsters, k = 2)
augment.gt_pca <- function(x, data = NULL, k = NULL, ...) {
  if (any(is.null(k), (k > ncol(x$u)))) {
    k <- ncol(x$u)
  }
  pred <- as.data.frame(sweep(x$u, 2, x$d, "*"))[, 1:k]
  names(pred) <- paste0(".fittedPC", seq_len(ncol(pred)))
  ret <- if (!missing(data) && !is.null(data)) {
    # check that names of the two columns are in sync
    if (!all.equal(data$id, rownames(as.data.frame(x$u)))) {
      stop(paste(
        "the data id column does not correspond to the individuals",
        "in the pca object 'x'"
      ))
    }
    if (inherits(data, "grouped_df")) {
      group_vars <- group_vars(data)
      data <- data %>% dplyr::ungroup()
      data <- data %>%
        dplyr::mutate(.rownames = data$id) %>%
        tibble::add_column(pred)
      data <- data %>% dplyr::group_by(across(all_of(group_vars)))
    } else {
      data %>%
        dplyr::mutate(.rownames = data$id) %>%
        tibble::add_column(pred)
    }
  } else {
    tibble(.rownames = rownames(as.data.frame(x$u[, 1:k]))) %>%
      add_column(pred)
  }

  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(ret)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(ret) <- obj_class
  }
  ret
}


#' Augment the loci table with information from a gt_pca object
#'
#' Augment for `gt_pca` accepts a model object and a `gen_tibble` and adds
#' loadings for each locus to the loci table. Loadings for each component are
#' stored in a separate column, which is given name with the pattern
#' ".loadingPC1", ".loadingPC2", etc. If `data` is missing, then a tibble with
#' the loadings is returned.
#' @param x  A `gt_pca` object returned by one of the `gt_pca_*` functions.
#' @param data the `gen_tibble` used to run the PCA.
#' @param k the number of components to add
#' @param ... Not used. Needed to match generic signature only.
#' @return A [gen_tibble] with a loadings added to the loci tibble (accessible
#'   with [show_loci()]. If `data` is missing, a tibble of loadings.
#' @export
#' @name augment_loci_gt_pca
#' @seealso [gt_pca_autoSVD()] [gt_pca_tidiers]
#' @examples
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
#' # Create PCA
#' pca <- gt_pca_partialSVD(lobsters)
#'
#' # Augment the gen_tibble with the PCA scores
#' augment_loci(pca, data = lobsters)
#'
augment_loci.gt_pca <- function(x, data = NULL, k = NULL, ...) {
  if (any(is.null(k), (k > ncol(x$v)))) {
    k <- ncol(x$v)
  }
  loadings <- as.data.frame(x$v)[, 1:k]
  names(loadings) <- paste0(".loadingPC", seq_len(ncol(loadings)))
  ret <- if (!missing(data) && !is.null(data)) {
    # check that names of the two columns are in sync
    # @TODO reinstate this check once we have rownames in the pca object
    # for loadings
    if (!identical(loci_names(data), rownames(as.data.frame(x$v)))) {
      stop(paste(
        "the loci names in 'data' do not correspond to the loci in",
        "the pca object 'x'"
      ))
    }
    show_loci(data) <- show_loci(data) %>% tibble::add_column(loadings)
  } else {
    # @TODO fir this once we have loci names in the pca object
    tibble(.rownames = rownames(x$v)) %>%
      add_column(loadings)
  }
  ret
}


# a print method
#' @method print gt_pca
#' @export
print.gt_pca <- function(x, ...) {
  cat(" === PCA of gen_tibble object ===")
  cat("\nMethod: ")
  print(x$method)
  cat("\nCall ($call):")
  print(x$call)
  cat(
    "\nEigenvalues ($d):\n",
    round(utils::head(x$d, 6), 3),
    ifelse(length(x$d) > 6, "...\n", "\n")
  )
  cat(
    "\nPrincipal component scores ($u):\n matrix with",
    nrow(x$u),
    "rows (individuals) and",
    ncol(x$u),
    "columns (axes)",
    "\n"
  )
  cat(
    "\nLoadings (Principal axes) ($v):\n matrix with",
    nrow(x$v),
    "rows (SNPs) and",
    ncol(x$v),
    "columns (axes)",
    "\n"
  )
  cat("\n")
}
