#' Print a summary of a merge report
#'
#' This function creates a summary of the merge report generated by
#' [rbind_dry_run()]
#'
#' @param object a list generated by [rbind_dry_run()]
#' @param ... unused (necessary for compatibility with generic function)
#' @param ref_label the label for the reference dataset (defaults to
#'   "reference")
#' @param target_label the label for the target dataset (defaults to "target")
#' @returns NULL (prints a summary to the console)
#' @rdname summary_rbind_dry_run
#' @aliases summary_rbind_report
#' @method summary rbind_report
#' @export

summary.rbind_report <- function(
    object,
    ...,
    ref_label = "reference",
    target_label = "target") {
  cat("harmonising loci between two datasets\n")
  cat(
    "flip_strand = ",
    attr(object, "flip_strand"),
    " ; remove_ambiguous = ",
    attr(object, "remove_ambiguous"),
    "\n"
  )
  cat("-----------------------------\n")
  cat("dataset:", ref_label, "\n")
  cat(
    "number of SNPs:",
    nrow(object$ref),
    "reduced to",
    sum(!is.na(object$ref$new_id)),
    "\n"
  )
  cat(
    "(",
    sum(object$ref$ambiguous),
    "are ambiguous, of which",
    (sum(object$ref$ambiguous & is.na(object$ref$new_id))),
    " were removed)\n"
  )
  cat("-----------------------------\n")
  cat("dataset:", target_label, "\n")
  cat(
    "number of SNPs:",
    nrow(object$target),
    "reduced to",
    sum(!is.na(object$target$new_id)),
    "\n"
  )
  cat(
    "(",
    sum(object$target$to_flip),
    "were flipped to match the reference set)\n"
  )
  cat(
    "(",
    sum(object$target$ambiguous),
    "are ambiguous, of which",
    (sum(object$target$ambiguous & is.na(object$target$new_id))),
    "were removed)"
  )
}
