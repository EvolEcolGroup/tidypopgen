#' Combine method for gt_admix objects
#'
#' @param ... A list of `gt_admix` objects
#' @param match_attributes boolean, determining whether all attributes (id,
#'  group and algorithm) of the `gt_admix` objects to be combined must be an
#'  exact match (TRUE, the default), or whether non-matching attributes should
#'  be ignored (FALSE)
#' @return A `gt_admix` object with the combined data
#' @export
#' @examples
#' # run the example only if we have the package installed
#' if (requireNamespace("LEA", quietly = TRUE)) {
#'   example_gt <- load_example_gt("gen_tbl")
#'
#'   # Create a gt_admix object
#'   admix_obj <- example_gt %>% gt_snmf(k = 1:3, project = "force")
#'
#'   # Create a second gt_admix object
#'   admix_obj2 <- example_gt %>% gt_snmf(k = 2:4, project = "force")
#'
#'   # Combine the two gt_admix objects
#'   new_admix_obj <- c(admix_obj, admix_obj2)
#'   summary(new_admix_obj)
#' }
c.gt_admix <- function(..., match_attributes = TRUE) {
  admix_list <- list(...)
  # check that all the objects are of class gt_admix
  if (!all(sapply(admix_list, function(x) inherits(x, "gt_admix")))) {
    stop("All the objects must be of class gt_admix")
  }

  # use first object as reference
  ref <- admix_list[[1]]
  ref_ids <- ref$id

  # check all Q matrices have same number of individuals
  n_inds <- nrow(ref$Q[[1]])
  for (obj in admix_list) {
    q_rows <- unique(sapply(obj$Q, nrow))
    if (length(q_rows) != 1 || q_rows != n_inds) {
      stop("All Q matrices must have the same number of individuals (rows)")
    }
  }

  # if all gt_admix objects have the same id,  check they are in the same order
  has_id <- sapply(admix_list, function(x) !is.null(x$id))

  if (all(has_id)) {
    ids_match <- all(sapply(admix_list, function(x) setequal(x$id, ref_ids)))
    order_match <- all(sapply(admix_list, function(x) identical(x$id, ref_ids)))

    if (ids_match && !order_match) {
      stop(paste(
        "All 'gt_admix' objects have the same individual ids, but individuals",
        "are in a different order. Please re-order individuals before",
        "combining objects using gt_admix_reorder_q()"
      ))
    }
  }

  # if match_attributes is TRUE, check that all attributes match
  if (match_attributes) {
    check_equal <- function(field) {
      all(sapply(admix_list, function(x) identical(x[[field]], ref[[field]])))
    }
    if (!check_equal("id")) stop("Id information does not match")
    if (!check_equal("group")) stop("Group information does not match")
    if (!check_equal("algorithm")) stop("Algorithm information does not match")
  }

  combined_obj <- list()

  # combine all the elements from each list
  combined_obj$k <- unlist(lapply(admix_list, function(x) x$k))
  combined_obj$Q <- do.call(c, lapply(admix_list, function(x) x$Q))

  # if we have a P element in all of the objects, combine it
  if (all(sapply(admix_list, function(x) !is.null(x$P)))) {
    combined_obj$P <- do.call(c, lapply(admix_list, function(x) x$P))
  }

  # if we have a log element in all of the objects, combine it
  if (all(sapply(admix_list, function(x) !is.null(x$log)))) {
    combined_obj$log <- do.call(c, lapply(admix_list, function(x) x$log))
  }

  # if we have a cv element in all of the objects, combine it
  if (all(sapply(admix_list, function(x) !is.null(x$cv)))) {
    combined_obj$cv <- unlist(lapply(admix_list, function(x) x$cv))
  }

  # if we have a log_lik element in any of the objects, combine it
  if (all(sapply(admix_list, function(x) !is.null(x$loglik)))) {
    combined_obj$loglik <- unlist(lapply(admix_list, function(x) x$loglik))
  }

  # handle id depending on match_attributes
  if (match_attributes) {
    # if match_attributes is TRUE, then ids are either all identical or all NULL
    if (!is.null(admix_list[[1]]$id)) {
      combined_obj$id <- admix_list[[1]]$id
    }
  } else {
    # if match_attributes is FALSE handle varied id presence and consistency
    id_list <- lapply(admix_list, function(x) x$id)
    has_id <- sapply(id_list, Negate(is.null))

    if (!any(has_id)) {
      # no IDs present at all, omit
    } else if (all(has_id) &&
      all(sapply(
        id_list[-1],
        function(x) identical(x, id_list[[1]])
      ))) {
      # all have id, and they match, use id from the first object
      combined_obj$id <- id_list[[1]]
    } else if (any(has_id)) {
      present_ids <- id_list[has_id]
      if (all(sapply(
        present_ids[-1],
        function(x) identical(x, present_ids[[1]])
      ))) {
        # some have id, all present ones match, use with warning
        combined_obj$id <- present_ids[[1]]
        warning("Only some gt_admix objects have id information; using this id")
      } else {
        # any mismatch among present ids, omit with warning
        warning("Id information is present but does not match. Omitting id")
      }
    }
  }

  # handle group depending on match_attributes
  if (match_attributes) {
    # if match_attributes is TRUE, then groups are either identical or all NULL
    if (!is.null(admix_list[[1]]$group)) {
      combined_obj$group <- admix_list[[1]]$group
    }
  } else {
    # if match_attributes is FALSE handle varied group presence and consistency
    group_list <- lapply(admix_list, function(x) x$group)
    has_group <- sapply(group_list, Negate(is.null))

    if (!any(has_group)) {
      # no groups present at all, omit
    } else if (all(has_group) &&
      all(sapply(
        group_list[-1],
        function(x) identical(x, group_list[[1]])
      ))) {
      # all have group, and they match, use group from the first object
      combined_obj$group <- group_list[[1]]
    } else if (any(has_group)) {
      present_groups <- group_list[has_group]
      if (all(sapply(
        present_groups[-1],
        function(x) identical(x, present_groups[[1]])
      ))) {
        # some have group, all present ones match, use with warning
        combined_obj$group <- present_groups[[1]]
        warning(paste(
          "Only some gt_admix objects have group information;",
          "using this group"
        ))
      } else {
        # any mismatch among present groups, omit with warning
        warning(paste(
          "Group information is present but does not match.",
          "Omitting group"
        ))
      }
    }
  }

  # handle algorithm depending on match_attributes
  if (match_attributes) {
    # if match_attributes is TRUE then algorithm is either identical or all NULL
    if (!is.null(admix_list[[1]]$algorithm)) {
      combined_obj$algorithm <- admix_list[[1]]$algorithm
    }
  } else {
    # if match_attributes=FALSE handle varied algorithm presence and consistency
    algorithm_list <- lapply(admix_list, function(x) x$algorithm)
    has_algorithm <- sapply(algorithm_list, Negate(is.null))

    if (!any(has_algorithm)) {
      # no algorithms present at all, omit
    } else if (all(has_algorithm) &&
      all(sapply(
        algorithm_list[-1],
        function(x) identical(x, algorithm_list[[1]])
      ))) {
      # all have algorithm, and they match, use algorithm from the first object
      combined_obj$algorithm <- algorithm_list[[1]]
    } else if (any(has_algorithm)) {
      present_algorithms <- algorithm_list[has_algorithm]
      if (all(sapply(
        present_algorithms[-1],
        function(x) identical(x, present_algorithms[[1]])
      ))) {
        # some have algorithm, all present ones match, use with warning
        combined_obj$algorithm <- present_algorithms[[1]]
        warning(paste(
          "Only some gt_admix objects have algorithm information;",
          "using this algorithm"
        ))
      } else {
        # any mismatch among algorithms, omit with warning
        warning(paste(
          "Algorithm information is present but does not match.",
          "Omitting algorithm"
        ))
      }
    }
  }

  # reorder combined object by k
  # get the order of k
  order_index <- order(combined_obj$k)

  # reorder all list elements by the new order
  combined_obj$Q <- combined_obj$Q[order_index]
  if (!is.null(combined_obj$P)) {
    combined_obj$P <- combined_obj$P[order_index]
  }
  if (!is.null(combined_obj$log)) {
    combined_obj$log <- combined_obj$log[order_index]
  }
  if (!is.null(combined_obj$cv)) {
    combined_obj$cv <- combined_obj$cv[order_index]
  }
  if (!is.null(combined_obj$loglik)) {
    combined_obj$loglik <- combined_obj$loglik[order_index]
  }

  # reorder k vector to match
  combined_obj$k <- combined_obj$k[order_index]

  # set the class of the object
  class(combined_obj) <- c("gt_admix", "list")
  return(combined_obj)
}

#' Summary method for gt_admix objects
#'
#' @param object a `gt_admix` object
#' @param ... unused (necessary for compatibility with generic function)
#' @return A summary of the `gt_admix` object
#' @export
#' @examples
#' # run the example only if we have the package installed
#' if (requireNamespace("LEA", quietly = TRUE)) {
#'   example_gt <- load_example_gt("gen_tbl")
#'
#'   # Create a gt_admix object
#'   admix_obj <- example_gt %>% gt_snmf(k = 1:3, project = "force")
#'
#'   # Print a summary
#'   summary(admix_obj)
#' }
#'
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
#' @examples
#' # run the example only if we have the package installed
#' if (requireNamespace("LEA", quietly = TRUE)) {
#'   example_gt <- load_example_gt("gen_tbl")
#'
#'   # Create a gt_admix object
#'   admix_obj <- example_gt %>% gt_snmf(k = 1:3, project = "force")
#'
#'   # The $id in admix_obj is the same as in the gen_tibble
#'   admix_obj$id
#'
#'   # Reorder the q matrices based on the grouping variable
#'   admix_obj <- gt_admix_reorder_q(admix_obj,
#'     group = example_gt$population
#'   )
#'
#'   # The $id in admix_obj is now reordered according to the population
#'   admix_obj$id
#' }
#'
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
