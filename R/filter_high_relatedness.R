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
        rownames(matrix) <- seq_len(nrow(matrix))
      } else {
        colnames(matrix) <- (.x)$id
        rownames(matrix) <- (.x)$id
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

    # take absolute value of each relationship
    matrix <- abs(matrix)

    # re-ordered columns based on max absolute correlation
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

    # loop through each individual
    # for pairs with higher than kings_threshold relatedness
    # the individual with higher average relatedness is removed
    for (i in 1:(var_num - 1)) {
      if (!any(matrix2[!is.na(matrix2)] > kings_threshold)) {
        if (verbose) {
          message("All correlations <=", kings_threshold, "\n")
        }
        break()
      }
      if (col_to_delete[i]) {
        next
      }
      for (j in (i + 1):var_num) {
        if (!col_to_delete[i] && !col_to_delete[j]) {
          if (matrix[i, j] > kings_threshold) {
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
              matrix2[i, ] <- NA
              matrix2[, i] <- NA
              if (verbose) {
                message(" so flagging column", new_order[i], "\n")
              }
            } else {
              col_to_delete[j] <- TRUE
              matrix2[j, ] <- NA
              matrix2[, j] <- NA
              if (verbose) {
                message(" so flagging column", new_order[j], "\n")
              }
            }
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
