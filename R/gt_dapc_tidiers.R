#' Tidy a `gt_dapc` object
#'
#' This summarizes information about the components of a `gt_dapc` from the
#' `tidypopgen` package. The parameter `matrix` determines which element is
#' returned.
#'
#' @param x A `gt_dapc` object (as returned by [gt_dapc()]).
#' @param matrix Character specifying which component of the DAPC should be
#'   tidied.
#'
#'   - `"samples"`, `"scores"`, or `"x"`: returns information about
#'   the map from the original space into the least discriminant axes.
#'
#'   - `"v"`, `"rotation"`, `"loadings"` or `"variables"`: returns information
#'   about the map from discriminant axes space back into the original space
#'   (i.e. the genotype frequencies). Note that this are different from the
#'   loadings linking to the PCA scores (which are available in the element
#'   $loadings of the dapc object).
#'
#'   - `"d"`, `"eigenvalues"` or `"lds"`: returns information about the
#'   eigenvalues.
#'
#' @param ... Not used. Needed to match generic signature only.
#'
#' @return A [tibble::tibble] with columns depending on the component of DAPC
#'   being tidied.
#'
#'   If `"scores"` each row in the tidied output corresponds to the original
#'   data in PCA space. The columns are:
#'
#'   \item{`row`}{ID of the original observation (i.e. rowname from original
#'     data).}
#'   \item{`LD`}{Integer indicating a principal component.}
#'   \item{`value`}{The score of the observation for that particular principal
#'     component. That is, the location of the observation in PCA space.}
#'
#'   If `matrix` is  `"loadings"`, each row in the tidied output corresponds to
#'   information about the principle components in the original space. The
#'   columns are:
#'
#'   \item{`row`}{The variable labels (colnames) of the data set on
#'   which PCA was performed.}
#'   \item{`LD`}{An integer vector indicating the principal component.}
#'   \item{`value`}{The value of the eigenvector (axis score) on the
#'   indicated principal component.}
#'
#'   If `"eigenvalues"`, the columns are:
#'
#'   \item{`LD`}{An integer vector indicating the discriminant axis.}
#'   \item{`std.dev`}{Standard deviation (i.e. sqrt(eig/(n-1))) explained by
#'   this DA (for compatibility with `prcomp`.}
#'   \item{`cumulative`}{Cumulative variation explained by
#'     principal components up to this component (note that this is NOT
#'     phrased as a percentage of total variance, since many methods
#'     only estimate a truncated SVD.}
#'
#' @aliases gt_dapc_tidiers
#' @export
#' @seealso [gt_dapc()] [augment.gt_dapc()]
#' @examples
#' #' # Create a gen_tibble of lobster genotypes
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
#' # Create PCA and run DAPC
#' pca <- gt_pca_partialSVD(lobsters)
#' populations <- as.factor(lobsters$population)
#' dapc_res <- gt_dapc(pca, n_pca = 6, n_da = 2, pop = populations)
#'
#' # Tidy scores
#' tidy(dapc_res, matrix = "scores")
#'
#' # Tidy eigenvalues
#' tidy(dapc_res, matrix = "eigenvalues")
#'
#' # Tidy loadings
#' tidy(dapc_res, matrix = "loadings")
#'
tidy.gt_dapc <- function(x, matrix = "eigenvalues", ...) {
  if (length(matrix) > 1) {
    stop("Must select a single matrix to tidy.", call. = FALSE)
  }

  matrix_args <- c(
    "rotation",
    "x",
    "variables",
    "samples",
    "v",
    "lds",
    "d",
    "scores",
    "loadings",
    "eigenvalues"
  )
  matrix <- rlang::arg_match(matrix, matrix_args)

  if (matrix %in% c("lds", "d", "eigenvalues")) {
    # square_frobenius <- sum(x$eig) #nolint
    ret <- tibble(
      LD = seq_len(length(x$eig)),
      "eigenvalue" = x$eig
    ) %>%
      # mutate(percent =  (std.dev)^2/nrow(x$d),
      mutate(cumulative = cumsum(.data$eigenvalue))
  } else if (matrix %in% c("rotation", "variables", "v", "loadings")) {
    ret <- x$var.load %>%
      tibble::as_tibble(rownames = "column") %>%
      tidyr::pivot_longer(
        cols = -"column",
        names_to = "LD",
        values_to = "value"
      )
  } else if (matrix %in% c("x", "samples", "scores")) {
    ret <- x$ind.coord
    ret <- ret %>%
      tibble::as_tibble(rownames = "row") %>%
      tidyr::pivot_longer(
        cols = -"row",
        names_to = "LD",
        values_to = "value"
      )
    # change the LD to a numeric
    ret <- mutate(
      ret,
      LD = as.numeric(str_replace_base(.data$LD, "LD", ""))
    )
  }

  as_tibble(ret)
}


#' Augment data with information from a gt_dapc object
#'
#' Augment for `gt_dapc` accepts a model object and a dataset and adds scores to
#' each observation in the dataset. Scores for each component are stored in a
#' separate column, which is given name with the pattern ".fittedLD1",
#' ".fittedLD2", etc. For consistency with [broom::augment.prcomp], a column
#' ".rownames" is also returned; it is a copy of 'id', but it ensures that any
#' scripts written for data augmented with [broom::augment.prcomp] will work out
#' of the box (this is especially helpful when adapting plotting scripts).
#' @param x  A `gt_dapc` object returned by [gt_dapc()].
#' @param data the `gen_tibble` used to run the PCA.
#' @param k the number of components to add
#' @param ... Not used. Needed to match generic signature only.
#' @return A  [gen_tibble] containing the original data along with additional
#'   columns containing each observation's projection into PCA space.
#' @export
#' @seealso [gt_dapc()] [gt_dapc_tidiers]
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
#' # Create PCA and run DAPC
#' pca <- gt_pca_partialSVD(lobsters)
#' populations <- as.factor(lobsters$population)
#' dapc_res <- gt_dapc(pca, n_pca = 6, n_da = 2, pop = populations)
#'
#' # Augment the gen_tibble with the DAPC scores
#' augment(dapc_res, data = lobsters, k = 2)
#'
augment.gt_dapc <- function(x, data = NULL, k = NULL, ...) {
  if (is.null(k)) {
    k <- ncol(x$ind.coord)
  }
  pred <- as.data.frame(x$ind.coord)[, 1:k]
  names(pred) <- paste0(".fittedLD", seq_len(ncol(pred)))
  ret <- if (!missing(data) && !is.null(data)) {
    # check that names of the two columns are in sync
    if (!all.equal(data$id, rownames(pred))) {
      stop(paste(
        "the data id column does not correspond to the individuals",
        "in the gt_dapc object 'x'"
      ))
    }
    data %>%
      dplyr::mutate(.rownames = data$id) %>%
      tibble::add_column(pred)
  } else {
    tibble(.rownames = rownames(pred)) %>%
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


# a print method
#' @method print gt_dapc
#' @export
print.gt_dapc <- function(x, ...) {
  cat(" === DAPC of gen_tibble object ===")
  cat("\nCall ($call):")
  print(x$call)
  cat(
    "\nEigenvalues ($eig):\n",
    round(utils::head(x$eig, 6), 3),
    ifelse(length(x$eig) > 6, "...\n", "\n")
  )
  cat(
    "\nLD scores ($ind.coord):\n matrix with",
    nrow(x$ind.coord),
    "rows (individuals) and",
    ncol(x$ind.coord),
    "columns (LD axes)",
    "\n"
  )
  cat(
    "\nLoadings by PC ($loadings):\n matrix with",
    nrow(x$loadings),
    "rows (PC axes) and",
    ncol(x$loadings),
    "columns (LD axes)",
    "\n"
  )
  if (!is.null(x$var.contr)) {
    cat(
      "\nLoadings by locus($var.load):\n matrix with",
      nrow(x$var.load),
      "rows (loci) and",
      ncol(x$var.load),
      "columns (LD axes)",
      "\n"
    )
  }
  cat("\n")
}

#
