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
  if (all(sapply(list(...), function(x) !is.null(x$log_lik)))) {
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
  return(combined_obj)
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
    tab_sum <- rbind(as.numeric(names(tab_sum)),tab_sum)
    rownames(tab_sum) <- c("k", "n")
    colnames(tab_sum) <- rep("", ncol(tab_sum))
    cat(" for multiple runs:")
    print(tab_sum)
  }
  cat("with slots:\n")
  cat("$Q for Q matrices\n")
  # if there is a lot P in the object, print it
  if ("P" %in% names(object)){
    cat("$P for  matrices\n")
  }
  # if there is a slot log in the object, print it
  if ("log" %in% names(object)){
    cat("$log for logs from the algorithm\n")
  }
  # if there is a slot cv in the object, print it
  if ("cv" %in% names(object)){
    cat("$cv for cross validation error\n")
  }
}
