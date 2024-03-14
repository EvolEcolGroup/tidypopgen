#' tidy a `gt_pca` object
#'
#' This summarizes information about the components of a `gt_pca` from the
#' `tidypopgen` package. The parameter `matrix` determines which element is
#' returned. Column names of the tidied output match those returned by
#' [broom::tidy.prcomp], the tidier for the standard pca objects returned
#' by [stats::prcomp].
#'
#' @param x A `gt_pca` object returned by [gt_pca].
#' @param matrix Character specifying which component of the PCA should be
#'   tidied.
#'
#'   - `"u"`, `"samples"`, `"scores"`, or `"x"`: returns information about
#'     the map from the original space into principle components space.
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
#'   \item{`std.dev`}{Standard deviation explained by this PC.}
#'   \item{`percent`}{Fraction of variation explained by this component
#'     (a numeric value between 0 and 1).}
#'   \item{`cumulative`}{Cumulative fraction of variation explained by
#'     principle components up to this component (a numeric value between 0 and
#'     1).}
#'
#' @aliases gt_pca_tidiers
#' @export
#' @seealso [gt_pca()] [augment.gt_pca()]

tidy.gt_pca <- function(x, matrix = "u", ...) {
  if (length(matrix) > 1) {
    stop("Must select a single matrix to tidy.", call. = FALSE)
  }

  MATRIX <- c(
    "rotation", "x", "variables", "samples", "v", "u", "pcs", "d",
    "scores", "loadings", "eigenvalues"
  )
  matrix <- rlang::arg_match(matrix, MATRIX)


  # totalvar = sum(scale^2)

  #ncomp <- NCOL(x$rotation)
  if (matrix %in% c("pcs", "d", "eigenvalues")) {
    total_var <- sum(sum(x$scale^2))
    ret <- tibble(PC = seq_len(length(x$d)),
                  "std.dev" = sqrt(x$d)) %>%
      mutate(percent = .data$std.dev/total_var,
             cumulative = cumsum(.data$percent))
  } else if (matrix %in% c("rotation", "variables", "v", "loadings")) {
    ret <- x$v %>%
      tibble::as_tibble(rownames = "column") %>%
      tidyr::pivot_longer(
        cols = -"column",
        names_to = "PC",
        values_to = "value"
      )
    if (is.null(rownames(x$rotation))) ret$column <- as.integer(ret$column)
  } else if (matrix %in% c("x", "samples", "u", "scores")) {
    ret <- x$u %>%
      tibble::as_tibble(rownames = "row") %>%
      tidyr::pivot_longer(
        cols = -"row",
        names_to = "PC",
        values_to = "value")
    if (is.null(rownames(x$scores))) ret$row <- as.integer(ret$row)
  }

  # change the PC to a numeric
  ret <- mutate(ret, PC = as.numeric(stringr::str_replace(.data$PC, "PC", "")))
  as_tibble(ret)
}


#' Augment data with information from a gt_pca object
#'
#' Augment for `gt_pca` accepts a model object and a dataset and adds
#' scores to each
#' observation in the dataset. Scores for each component are stored in a
#' separate column, which is given name with the pattern ".fittedPC1",
#' ".fittedPC2", etc. For consistency with [broom::augment.prcomp], a column
#' ".rownames" is also returned; it is a copy of 'id', but it ensures that
#' any scripts written for data augmented with [broom::augment.prcomp] will
#' work out of the box (this is especially helpful when adapting plotting scripts).
#' @param x A `gt_pca` object returned by [gt_pca()].
#' @param data the `gen_tibble` used to run [gt_pca()].
#' @param k the number of components to add
#' @param ... Not used. Needed to match generic signature only.
#' @return A  [gen_tibble] containing the original data along with
#'   additional columns containing each observation's projection into
#'   PCA space.
#' @export
#' @seealso [gt_pca()] [gt_pca_tidiers]

augment.gt_pca <- function(x, data = NULL, k= NULL, ...) {
    if (is.null(k) | (k > ncol(x$u))){
      k <- ncol(x$u)
    }
    pred <- as.data.frame(x$u[,1:k])
    names(pred) <- paste0(".fittedPC", seq_len(ncol(pred)))
    # browser()
    ret <- if (!missing(data) && !is.null(data)) {
      #check that names of the two columns are in sync
      if (!all.equal(data$id, rownames(as.data.frame(x$u)))){
        stop("the data id column does not correspond to the individuals in the pca object 'x'")
      }
      data %>% tibble::add_column(.rownames = data$id, pred)
    } else {
      tibble(.rownames = rownames(as.data.frame(x$u[,1:k]))) %>%
        add_column(pred)
    }
    ret
}


