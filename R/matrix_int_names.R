#' Matrix with integer row/column names class
#'
#' Create a matrix with integer or character names for rows and columns. Integer
#' names are useful when dealing with very large numbers of rows or columns, as
#' the integer names can be stored more efficiently than character names. Whilst
#' this is possible for `data.frames`, the native `matrix` class only allows for
#' character row and column names via the `dimnames` attribute. Use
#' [row_names()] and [col_names()] to get or set the integer names, which are
#' stored in special attributes (`int_rownames` and `int_colnames`) on the
#' matrix.
#'
#' This class allows you to have integer names as attributes while still being a
#' matrix, and provides methods for getting and setting these names. When you
#' create a `matrix_int_names` object, you can provide integer or character
#' vectors for row and column names. If you provide integer vectors, they will
#' be stored in special attributes (`int_rownames` and `int_colnames`) instead
#' of the standard `dimnames`. If you provide character vectors, they will be
#' stored in the standard `dimnames` as usual. You can also mix and match,
#' having integer names for rows and character names for columns, or vice versa.
#' Note that, since the row and column names are stored in special attributes,
#' you have to use [row_names()] and [col_names()] to get or set them, rather
#' than `rownames()` and `colnames()`, which will only return character names if
#' present.
#'
#' @param data Matrix data or object coercible to matrix
#' @param row_names Integer vector (stored as int_rownames) or character vector
#'   (stored as dimnames)
#' @param col_names Integer vector (stored as int_colnames) or character vector
#'   (stored as dimnames)
#' @return A `matrix_int_names` object
#' @family matrix_int_names_functions
#' @export
#' @examples
#' # Create a matrix with integer row and column names
#' my_mat <- matrix_int_names(matrix(1:6, nrow = 2),
#'  row_names = c(10L, 20L),
#'  col_names = c(100L, 200L, 300L)
#' )
#' row_names(my_mat) # returns integer row names
#' col_names(my_mat) # returns integer column names
matrix_int_names <- function(data, row_names = NULL, col_names = NULL) {
  # Convert to matrix if needed
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  # Set class
  class(data) <- c("matrix_int_names", "matrix", "array")

  # Handle row names
  if (!is.null(row_names)) {
    row_names(data) <- row_names
  }

  # Handle column names
  if (!is.null(col_names)) {
    col_names(data) <- col_names
  }

  return(data)
}


#' Get or set column names for matrix_int_names
#'
#' This function reliably returns integer names if present, otherwise character
#' dimnames. Use this instead of colnames() if you need guaranteed dispatch.
#'
#' @param x A matrix_int_names object
#' @return Integer or character vector of column names, or NULL
#' @family matrix_int_names_functions
#' @export
col_names <- function(x) {
  UseMethod("col_names")
}

#' @export
#' @rdname col_names
col_names.default <- function(x) {
  base::colnames(x)
}

#' Get column names (both integer and character)
#' @export
#' @rdname col_names
col_names.matrix_int_names <- function(x) {
  # if we have the attribute for integer column names, return that, otherwise
  # return character column names
  if (!is.null(attr(x, "int_colnames"))) {
    return(attr(x, "int_colnames"))
  } else {
    return(colnames(x))
  }
}

#' @param value a vector of integer or character column names to set, replacing
#'   the current one.
#' @export
#' @rdname col_names
"col_names<-" <- function(x, value) {
  UseMethod("col_names<-", x)
}

#' default method
#' @export
#' @rdname col_names
`col_names<-.default` <- function(x, value) {
  base::colnames(x) <- value
  x
}

#' Set column names (integer or character)
#' @export
#' @rdname col_names
`col_names<-.matrix_int_names` <- function(x, value) {
  if (is.null(value)) {
    # Clear both types of names
    attr(x, "int_colnames") <- NULL
    class(x) <- c("matrix", "array")
    colnames(x) <- NULL
    class(x) <- c("matrix_int_names", "matrix", "array")
  } else if (is.numeric(value)) {
    if (!is.integer(value)) {
      # test if these are integer-valued, and thus can be coerced to integer
      if (is_integer_valued(value)) {
        value <- as.integer(round(value))
      } else {
        stop("col_names must be integer or character vector")
      }
    }    # Set as integer names
    if (length(value) != ncol(x)) {
      stop("Length of col_names must match number of columns")
    }
    attr(x, "int_colnames") <- value
    # Clear character dimnames for columns
    dn <- dimnames(x)
    if (!is.null(dn)) {
      dn[[2]] <- NULL
      dimnames(x) <- dn
    }
  } else if (is.character(value)) {
    # Set as character names
    if (length(value) != ncol(x)) {
      stop("Length of col_names must match number of columns")
    }
    # Clear integer colnames
    attr(x, "int_colnames") <- NULL
    # Set character dimnames
    colnames(x) <- value
  } else {
    stop("col_names must be integer or character vector")
  }
  x
}

#' Get or set row names for matrix_int_names
#'
#' Returns integer names if present, otherwise character dimnames. Use instead
#' of `rownames()` for guaranteed dispatch.
#'
#' @param x A matrix_int_names object
#' @family matrix_int_names_functions
#' @return Integer or character vector of row names, or NULL

#' @export
row_names <- function(x) {
  UseMethod("row_names")
}

#' @export
row_names.default <- function(x) {
  base::rownames(x)
}


#'
#' @export
row_names.matrix_int_names <- function(x) {
  int_rnames <- attr(x, "int_rownames")

  # Return integer names if present
  if (!is.null(int_rnames)) {
    return(int_rnames)
  } else {
    # Otherwise get character dimnames
    rownames(x)
  }
}

#' Set row names for matrix_int_names
#'
#' @param x A matrix_int_names object
#' @param value Integer or character vector of row names, or NULL
#' @return The modified matrix_int_names object
#' @export
#' @family matrix_int_names_functions
#' @rdname row_names
"row_names<-" <- function(x, value) {
  UseMethod("row_names<-", x)
}

# default method
#' @rdname row_names
#' @export
`row_names<-.default` <- function(x, value) {
  base::rownames(x) <- value
  x
}

#' @rdname row_names
#' @export
`row_names<-.matrix_int_names` <- function(x, value) {
  if (is.null(value)) {
    # Clear both types of names
    attr(x, "int_rownames") <- NULL
    dn <- dimnames(x)
    if (!is.null(dn)) {
      dn[[1]] <- NULL
      dimnames(x) <- dn
    }
  } else if (is.numeric(value)) {
    if (!is.integer(value)) {
      # test if these are integer-valued, and thus can be coerced to integer
      if (is_integer_valued(value)) {
        value <- as.integer(round(value))
      } else {
        stop("row_names must be integer or character vector")
      }
    }

    # Set as integer names
    if (length(value) != nrow(x)) {
      stop("Length of row_names must match number of rows")
    }
    attr(x, "int_rownames") <- value
    # Clear character dimnames for rows
    dn <- dimnames(x)
    if (!is.null(dn)) {
      dn[[1]] <- NULL
      dimnames(x) <- dn
    }
  } else if (is.character(value)) {
    # Set as character names
    if (length(value) != nrow(x)) {
      stop("Length of row_names must match number of rows")
    }
    # Clear integer rownames
    attr(x, "int_rownames") <- NULL
    # Set character dimnames
    dn <- dimnames(x)
    if (is.null(dn)) {
      dn <- list(NULL, NULL)
    }
    dn[[1]] <- value
    dimnames(x) <- dn
  } else {
    stop("row_names must be integer or character vector")
  }
  x
}

# a small internal function to check if a vector is integer-valued (i.e. all
# values are close to integers within a tolerance)
is_integer_valued <- function(x, tol = .Machine$double.eps^0.5) {
  is.numeric(x) && all(abs(x - round(x)) < tol, na.rm = TRUE)
}



#' Subsetting method for matrix_int_names
#' @param x A matrix_int_names object
#' @param i Row indices (numeric), row names (character), or missing)
#' @param j Column indices (numeric), column names (character), or missing)
#' @param i_names Optional integer row names to subset by
#' @param j_names Optional integer column names to subset by
#' @param drop Logical indicating whether to drop dimensions
#' @return A subsetted matrix_int_names object
#' @export
#' @family matrix_int_names_functions
#' @rdname subset
`[.matrix_int_names` <- function(x, i, j, i_names = NULL,
                                 j_names = NULL, drop = TRUE) {
  # check that we don't have i and i_names both provided
  if (!missing(i) && !is.null(i_names)) {
    stop("Provide either i or i_names, not both")
  }
  # same for j
  if (!missing(j) && !is.null(j_names)) {
    stop("Provide either j or j_names, not both")
  }

  # Get the integer names and character dimnames
  rnames <- attr(x, "int_rownames")
  cnames <- attr(x, "int_colnames")

  # Get character dimnames directly
  char_rnames <- dimnames(x)[[1]]
  char_cnames <- dimnames(x)[[2]]

  # if i is missing, set to all rows or use i_names
  if (missing(i)) {
    # if we have i_names, use them
    if (!is.null(i_names)) {
      # if i_names are characters, convert to indices
      if (is.character(i_names)) {
        i <- i_names # we will convert later
      } else {
        # Check if these are position indices or name indices
        if (all(i_names %in% rnames)) {
          i <- match(i_names, rnames)
        } else {
          stop("some i_names do not match any integer row names")
        }
      }
    } else {
      i <- seq_len(nrow(x))
    }
  }

  # do the same for j
  if (missing(j)) {
    # if we have j_names, use them
    if (!is.null(j_names)) {
      if (is.character(j_names)) {
        j <- j_names # we will convert later
      } else {
        # Check if these are position indices or name indices
        if (all(j_names %in% cnames)) {
          j <- match(j_names, cnames)
        } else {
          stop("some j_names do not match any integer column names")
        }
      }
    } else {
      j <- seq_len(ncol(x))
    }
  }


  # Convert character names to indices if needed
  if (!missing(i) && is.character(i) && !is.null(char_rnames)) {
    i <- match(i, char_rnames)
  }

  if (!missing(j) && is.character(j) && !is.null(char_cnames)) {
    j <- match(j, char_cnames)
  }

  class(x) <- c("matrix", "array") # Temporarily remove class
  # Perform subsetting on the underlying matrix
  result <- x[i, j, drop = drop]

  # Preserve matrix_int_names class and update attributes if not dropped to
  # vector
  if (is.matrix(result)) {
    # Update integer names if present
    new_rnames <- if (!is.null(rnames)) rnames[i] else NULL
    new_cnames <- if (!is.null(cnames)) cnames[j] else NULL

    # Update character names if present
    new_char_rnames <- if (!is.null(char_rnames)) char_rnames[i] else NULL
    new_char_cnames <- if (!is.null(char_cnames)) char_cnames[j] else NULL

    # Set integer names
    if (!is.null(new_rnames)) {
      attr(result, "int_rownames") <- new_rnames
    }
    if (!is.null(new_cnames)) {
      attr(result, "int_colnames") <- new_cnames
    }

    # Set character dimnames
    if (!is.null(new_char_rnames) || !is.null(new_char_cnames)) {
      dimnames(result) <- list(new_char_rnames, new_char_cnames)
    }

    # Restore class
    class(result) <- c("matrix_int_names", "matrix", "array")
  }

  return(result)
}

#' Assignment method for matrix_int_names
#' @param value A value to assign to the subsetted matrix
#' @rdname subset
#' @export
`[<-.matrix_int_names` <- function(x, i, j, i_names = NULL,
                                   j_names = NULL, value) {
  # Get the integer names and character dimnames
  rnames <- attr(x, "int_rownames")
  cnames <- attr(x, "int_colnames")

  # Get character dimnames directly
  char_rnames <- dimnames(x)[[1]]
  char_cnames <- dimnames(x)[[2]]

  # if i is missing, set to all rows or use i_names
  if (missing(i)) {
    # if we have i_names, use them
    if (!is.null(i_names)) {
      # if i_names are characters, convert to indices
      if (is.character(i_names)) {
        i <- i_names # we will convert later
      } else {
        # Check if these are position indices or name indices
        if (all(i_names %in% rnames)) {
          i <- match(i_names, rnames)
        } else {
          stop("some i_names do not match any integer row names")
        }
      }
    } else {
      i <- seq_len(nrow(x))
    }
  }

  # do the same for j
  if (missing(j)) {
    # if we have j_names, use them
    if (!is.null(j_names)) {
      if (is.character(j_names)) {
        j <- j_names # we will convert later
      } else {
        # Check if these are position indices or name indices
        if (all(j_names %in% cnames)) {
          j <- match(j_names, cnames)
        } else {
          stop("some j_names do not match any integer column names")
        }
      }
    } else {
      j <- seq_len(ncol(x))
    }
  }

  # Convert character names to indices if needed
  if (!missing(i) && is.character(i) && !is.null(char_rnames)) {
    i <- match(i, char_rnames)
  }

  if (!missing(j) && is.character(j) && !is.null(char_cnames)) {
    j <- match(j, char_cnames)
  }

  # Perform assignment
  class(x) <- c("matrix", "array") # Temporarily remove matrix_int_names class
  x[i, j] <- value

  # Restore class and attributes
  class(x) <- c("matrix_int_names", "matrix", "array")

  # DEBUG this is probably not needed
  attr(x, "int_rownames") <- rnames
  attr(x, "int_colnames") <- cnames
  if (!is.null(char_rnames) || !is.null(char_cnames)) {
    dimnames(x) <- list(char_rnames, char_cnames)
  }

  return(x)
}

#' Print method for matrix_int_names
#'
#' @param x A matrix_int_names object
#' @param ... Additional arguments passed to print
#' @return Invisibly returns the original object
#' @family matrix_int_names_functions
#' @export
print.matrix_int_names <- function(x, ...) {
  cat("Matrix with integer/character names\n")
  cat("Dimensions:", nrow(x), "x", ncol(x), "\n")

  rnames <- attr(x, "int_rownames")
  cnames <- attr(x, "int_colnames")

  # Get character dimnames directly
  char_rnames <- dimnames(x)[[1]]
  char_cnames <- dimnames(x)[[2]]

  if (!is.null(rnames)) {
    cat("Integer row names:", paste(rnames, collapse = ", "), "\n")
  }
  if (!is.null(cnames)) {
    cat("Integer column names:", paste(cnames, collapse = ", "), "\n")
  }
  if (!is.null(char_rnames)) {
    cat("Character row names:", paste(char_rnames, collapse = ", "), "\n")
  }
  if (!is.null(char_cnames)) {
    cat("Character column names:", paste(char_cnames, collapse = ", "), "\n")
  }

  cat("\n")

  # Create a display version
  display <- as.matrix(x)
  class(display) <- c("matrix", "array")

  # Prefer integer names for display, fallback to character
  if (!is.null(rnames)) {
    dn <- dimnames(display)
    if (is.null(dn)) dn <- list(NULL, NULL)
    dn[[1]] <- as.character(rnames)
    dimnames(display) <- dn
  } else if (!is.null(char_rnames)) {
    dn <- dimnames(display)
    if (is.null(dn)) dn <- list(NULL, NULL)
    dn[[1]] <- char_rnames
    dimnames(display) <- dn
  }

  if (!is.null(cnames)) {
    dn <- dimnames(display)
    if (is.null(dn)) dn <- list(NULL, NULL)
    dn[[2]] <- as.character(cnames)
    dimnames(display) <- dn
  } else if (!is.null(char_cnames)) {
    dn <- dimnames(display)
    if (is.null(dn)) dn <- list(NULL, NULL)
    dn[[2]] <- char_cnames
    dimnames(display) <- dn
  }

  print(display, ...)
  invisible(x)
}


#' Coercion method to data.frame for matrix_int_names
#'
#' This method converts a matrix_int_names object to a data.frame, preserving
#' the integer row names. Otherwise, it works similarly to the default
#' as.data.frame.matrix method.
#' 
#' @param x A matrix_int_names object
#' @param row.names Ignored, as row names are preserved from the matrix_int_names
#' @param optional Ignored, included for compatibility with generic signature
#' @param ... Additional arguments passed to as.data.frame.matrix
#' @export
#' @family matrix_int_names_functions
#' @examples
#' # Create a matrix with integer row and column names
#' my_mat <- matrix_int_names(matrix(1:6, nrow = 2),
#'   row_names = c(10L, 20L),
#'   col_names = c(100L, 200L, 300L)
#' )
#' # Convert to data.frame
#' my_df <- as.data.frame(my_mat)
#' my_df
#' 
as.data.frame.matrix_int_names <- function(x, row.names = NULL,
                                           optional = FALSE,
                                           ...) {
  df <- as.data.frame.matrix(x, row.names = row.names, optional = optional, ...)
  if (is.null(row.names)){
    attr(df, "row.names") <- row_names(x)
  }
  names(df) <- col_names(x)
  return(df)
}