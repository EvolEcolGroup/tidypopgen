#' Combine method for gt_admix objects
#'
#' @param ... A list of `gt_admix` objects
#' @return A `gt_admix` object with the combined data
#' @export

c.gt_admix <- function(...) {
  # check that all the objects are of class gt_admix
  if (!all(sapply(list(...), function(x) inherits(x, "gt_admix")))) {
    stop("All the objects must be of class gt_admix")
  }

  combined_obj <- list()
  # combine all the elements from each list
  combined_obj$k <- sapply(list(...), function(x) x$k)
  combined_obj$Q <- sapply(list(...), function(x) x$Q)
  # if we have a P element in any of the objects, combine it
  if (all(sapply(list(...), function(x) !is.null(x$P)))) {
    combined_obj$P <- sapply(list(...), function(x) x$P)
  }
  # if we have a log_lik element in any of the objects, combine it
  if (all(sapply(list(...), function(x) !is.null(x$loglik)))) {
    combined_obj$loglik <- unlist(sapply(list(...), function(x) x$loglik))
  }

  # if we have a log element in any of the objects, combine it
  if (all(sapply(list(...), function(x) !is.null(x$log)))) {
    combined_obj$log <- sapply(list(...), function(x) x$log)
  }
  # if we have a cv element in any of the objects, combine it
  if (all(sapply(list(...), function(x) !is.null(x$cv)))) {
    combined_obj$cv <- unlist(sapply(list(...), function(x) x$cv))
  }
  # if the first object has an id element, use it in the combined object
  if (!is.null(list(...)[[1]]$id)) {
    combined_obj$id <- list(...)[[1]]$id
  }
  # if the first object has a group element, use it in the combined object
  if (!is.null(list(...)[[1]]$group)) {
    combined_obj$group <- list(...)[[1]]$group
  }
  # set the class of the object
  class(combined_obj) <- c("gt_admix", "list")
  return(combined_obj) # nolint
}

#' Summary method for gt_admix objects
#'
#' @param object a `gt_admix` object
#' @param ... unused (necessary for compatibility with generic function)
#' @return A summary of the `gt_admix` object
#' @export
summary.gt_admix <- function(object, ...) {
  cat("Admixture results")
  # if we only have one element, give the k
  if (length(object$k) == 1) {
    cat(" for k = ", object$k, "\n")
  } else {
    tab_sum <- table(object$k)
    tab_sum <- rbind(as.numeric(names(tab_sum)), tab_sum)
    rownames(tab_sum) <- c("k", "n")
    colnames(tab_sum) <- rep("", ncol(tab_sum))
    cat(" for multiple runs:")
    print(tab_sum)
  }
  cat("with slots:\n")
  cat("$Q for Q matrices\n")
  # if there is a lot P in the object, print it
  if ("P" %in% names(object)) {
    cat("$P for  matrices\n")
  }
  # if there is a slot log in the object, print it
  if ("log" %in% names(object)) {
    cat("$log for logs from the algorithm\n")
  }
  # if there is a slot cv in the object, print it
  if ("cv" %in% names(object)) {
    cat("$cv for cross validation error\n")
  }
}


#' Reorder the q matrices based on the grouping variable
#'
#' This function reorders the q matrices in a `gt_admix` object based on the
#' grouping variable. This is useful before plotting when the samples from each
#' group are not adjacent to each other in the q matrix.
#'
#' @param x a `gt_admix` object, possibly with a grouping variable
#' @param group a character vector with the grouping variable (if there is no
#'   grouping variable info in `x`)
#' @return a `gt_admix` object with the q matrices reordered
#' @export

gt_admix_reorder_q <- function(x, group = NULL) {
  # check that x is a gt_admix object
  if (!inherits(x, "gt_admix")) {
    stop("x must be a gt_admix object")
  }
  # if we have a gruop variable,
  if (!is.null(group)) {
    # check that it is the same length as the q matrix
    if (length(group) != nrow(x$Q[[1]])) {
      stop(paste(
        "The length of the group variable must be the same as the",
        "number of rows in the Q matrix"
      ))
    }
    # and use it
    if (!is.null(group)) {
      x$group <- group
    }
  }
  # if we have no group variable, check that we have one in the object
  if (is.null(x$group)) {
    # if group is null
    if (is.null(group)) {
      stop(paste(
        "You must provide a group variable if there is no grouping",
        "information in the gt_admix object"
      ))
    }
  }

  group_meta <- tibble(id = seq(1, length(x$group)), group = x$group)
  # sort group meta by group
  group_meta <- group_meta %>% arrange(.data$group)
  # reorder the q matrices
  x$Q <- lapply(x$Q, function(y) y[group_meta$id, ])
  # if we have an id element, reorder it
  if (!is.null(x$id)) {
    x$id <- x$id[group_meta$id]
  } else {
    x$id <- group_meta$id
  }
  # reorder the group element
  x$group <- x$group[group_meta$id]

  return(x)
}
