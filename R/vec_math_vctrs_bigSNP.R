#' Mathematical predicates for vctrs_bigSNP vectors
#'
#' Provides support for the mathematical predicate functions
#' [is.finite()], [is.infinite()], and [is.nan()] on
#' `vctrs_bigSNP` vectors.
#'
#' These methods are primarily required for compatibility with
#' RStudio's object inspector and other tooling that may call
#' mathematical predicates on custom vctrs classes.
#'
#' All other Math-group generics are deliberately unsupported and
#' will raise an error.
#'
#' @param .fn Name of the mathematical function being dispatched.
#' @param .x A `vctrs_bigSNP` vector.
#' @param ... Unused additional arguments.
#'
#' @return
#' A logical vector for `is.finite()`, `is.infinite()`, and
#' `is.nan()`.
#'
#' @export
vec_math.vctrs_bigSNP <- function(.fn, .x, ...) { # nolint

  # Only allow the three predicate functions that RStudio
  # commonly uses when inspecting objects.
  if (!.fn %in% c("is.finite", "is.infinite", "is.nan")) {
    stop(
      sprintf(
        "`%s()` is not implemented for <vctrs_bigSNP>.",
        .fn
      ),
      call. = FALSE
    )
  }

  # Extract the underlying numeric data.
  data <- vctrs::vec_data(.x)

  # Dispatch to the corresponding base R function.
  switch(.fn,
    "is.finite"   = is.finite(data),
    "is.infinite" = is.infinite(data),
    "is.nan"      = is.nan(data)
  )
}
