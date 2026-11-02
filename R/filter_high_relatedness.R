#' Filter individuals based on a relationship threshold
#'
#' This function takes a matrix of x by y individuals containing relatedness
#' coefficients and returns the maximum set of individuals that contains no
#' relationships above the given threshold.
#'
#'
#' @param matrix a square symmetric matrix of individuals containing
#'   relationship coefficients
#' @param .x a [`gen_tibble`] object
#' @param kings_threshold a threshold over which
#' @param verbose boolean whether to report to screen
#' @details Missing relationship coefficients (e.g. KING kinship for pairs of
#'   individuals with no loci genotyped in both) are not evidence of
#'   relatedness: such pairs are treated as not related, with a warning.
#' @return a list where '1' is individual ID's to retain, '2' is individual ID's
#'   to remove, and '3' is a boolean where individuals to keep are TRUE and
#'   individuals to remove are FALSE
#' @rdname filter_high_relatedness
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Calculate relationship matrix
#' king_matrix <- example_gt %>% pairwise_king(as_matrix = TRUE)
#'
#' # Filter individuals with threshold above 0.2
#' filter_high_relatedness(king_matrix, example_gt, kings_threshold = 0.2)
filter_high_relatedness <-
  function(matrix, .x = NULL, kings_threshold = NULL, verbose = FALSE) {
    # get number of individuals
    var_num <- dim(matrix)[1]

    # append row and col names
    if (is.null(dimnames(matrix))) {
      if (is.null(.x)) {
        colnames(matrix) <- seq_len(ncol(matrix))
        row.names(matrix) <- seq_len(nrow(matrix))
      } else {
        colnames(matrix) <- (.x)$id
        row.names(matrix) <- (.x)$id
      }
    }

    # get individual names
    var_names <- dimnames(matrix)[[1]]

    # if there is only 1 individual
    if (var_num == 1) {
      # return variable names
      passed_filter <- var_names
      to_remove <- character(0)
      var_names <- var_names %in% passed_filter == TRUE
      return(list(passed_filter, to_remove, var_names))
    }

    # missing kinship (e.g. no loci genotyped in both individuals) is not
    # evidence of relatedness, so such pairs are treated as not related; we
    # warn, as these pairs could not be compared
    n_missing <- sum(is.na(matrix[upper.tri(matrix)]))
    if (n_missing > 0) {
      warning(
        n_missing, " pairs of individuals have missing kinship (e.g. no loci ",
        "genotyped in both), so they could not be compared, and are not ",
        "treated as related"
      )
    }

    # note that we do not take the absolute value of each relationship, as
    # negative kinship coefficients indicate unrelated individuals

    # re-ordered columns based on average relatedness
    original_order <- 1:var_num

    # function to calculate average relatedness
    average_corr <-
      function(matrix) {
        mean(matrix, na.rm = TRUE)
      }
    tmp <- matrix

    # remove self-relatedness (diagonal)
    diag(tmp) <- NA

    # calculate average relatedness and a new order
    max_abs_cor_order <-
      order(apply(tmp, 2, average_corr), decreasing = TRUE)

    # re-order individuals in matrix based on relatedness
    matrix <- matrix[max_abs_cor_order, max_abs_cor_order]

    # record new order
    new_order <- original_order[max_abs_cor_order]
    rm(tmp)

    # initialize new variables
    col_to_delete <- rep(FALSE, var_num)
    matrix2 <- matrix
    diag(matrix2) <- NA

    # which relationships are above the threshold (missing kinship, e.g. no
    # loci genotyped in both individuals, is not evidence of relatedness)
    is_above <- function(x) !is.na(x) & x > kings_threshold
    # number of relationships above the threshold among the individuals that
    # have not been removed; rather than scanning the whole matrix at every
    # iteration (which costs O(n^3) overall), this count is updated whenever
    # an individual is removed
    n_above <- sum(is_above(matrix2))
    # number of relationships above the threshold in row and column k of
    # matrix2 (i.e. those lost when individual k is removed; the diagonal is
    # NA, so nothing is counted twice)
    n_above_k <- function(matrix2, k) {
      sum(is_above(matrix2[k, ])) + sum(is_above(matrix2[, k]))
    }

    # loop through each individual
    # for pairs with higher than kings_threshold relatedness
    # the individual with higher average relatedness is removed
    for (i in 1:(var_num - 1)) {
      if (n_above == 0) {
        if (verbose) {
          message("All correlations <=", kings_threshold, "\n")
        }
        break()
      }
      if (col_to_delete[i]) {
        next
      }
      # the individuals after i (that have not been removed) whose
      # relationship with i is above the threshold; removing one of them
      # does not affect the others, so they can all be found at once
      js <- (i + 1):var_num
      js <- js[!col_to_delete[js] & is_above(matrix[i, js])]
      for (j in js) {
        # stop comparing i once it has been removed
        if (col_to_delete[i]) {
          break()
        }
        mn1 <- mean(matrix2[i, ], na.rm = TRUE)
        mn2 <- mean(matrix2[-j, ], na.rm = TRUE)
        if (verbose) {
          message(
            "Compare row",
            new_order[i],
            " and column ",
            new_order[j],
            "with corr ",
            round(matrix[i, j], 3),
            "\n"
          )
        }
        if (verbose) {
          message("  Means: ", round(mn1, 3), "vs", round(mn2, 3))
        }
        if (mn1 > mn2) {
          col_to_delete[i] <- TRUE
          n_above <- n_above - n_above_k(matrix2, i)
          matrix2[i, ] <- NA
          matrix2[, i] <- NA
          if (verbose) {
            message(" so flagging column", new_order[i], "\n")
          }
        } else {
          col_to_delete[j] <- TRUE
          n_above <- n_above - n_above_k(matrix2, j)
          matrix2[j, ] <- NA
          matrix2[, j] <- NA
          if (verbose) {
            message(" so flagging column", new_order[j], "\n")
          }
        }
      }
    }

    # return variable names
    passed_filter <- var_names[new_order][!col_to_delete]
    to_remove <- var_names[!var_names %in% passed_filter]

    var_names <- var_names %in% passed_filter == TRUE

    return(list(passed_filter, to_remove, var_names))
  }
