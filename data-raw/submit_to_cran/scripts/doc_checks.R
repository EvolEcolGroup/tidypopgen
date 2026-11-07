# Checks of the Rd files of an R package that CRAN reviewers ask about, but
# that R CMD check does not report: every help page of a function needs a
# \value section (roxygen @returns) and an \examples section (roxygen
# @examples).
#
# Usage, from any working directory:
#   source("<skill directory>/scripts/doc_checks.R")
#   check_returns_in_documentation("path/to/package")
#   check_examples_in_documentation("path/to/package")
#
# Only base R is used (the Rd files are parsed with tools::parse_Rd(), rather
# than searched with grep), so the checks also work on Windows.

#' The top-level sections of an Rd file
#'
#' @param rd_file the path to an Rd file.
#' @returns a character vector with the tags of the top-level sections of the
#'   file (e.g. `"\\usage"`, `"\\value"`, `"\\examples"`), and, for
#'   `\docType{}` and `\keyword{}`, the tag followed by its value (e.g.
#'   `"\\docType{data}"`, `"\\keyword{internal}"`).
#' @keywords internal
rd_sections <- function(rd_file) {
  rd <- tools::parse_Rd(rd_file, permissive = TRUE)
  tags <- vapply(rd, function(x) attr(x, "Rd_tag"), character(1))
  with_value <- tags %in% c("\\docType", "\\keyword")
  values <- vapply(
    which(with_value),
    function(i) {
      paste0(tags[i], "{", trimws(paste(unlist(rd[[i]]), collapse = "")), "}")
    },
    character(1)
  )
  c(tags, values)
}

#' The Rd files that document functions
#'
#' Help pages that do not document functions do not need `\value` or
#' `\examples`, so they are left out:
#' - datasets (`\docType{data}`), and any other page without `\usage` (e.g.
#'   the documentation of external data);
#' - the package page (`\docType{package}`);
#' - re-exported functions (`\docType{import}`);
#' - if `skip_internal` is `TRUE`, pages with `\keyword{internal}` (roxygen
#'   `@keywords internal`), which are not listed in the index of the help;
#' - the pages given in `exclude`.
#'
#' @param pkg the path to the root of the package (the directory with the
#'   DESCRIPTION file).
#' @param exclude the names of Rd files (e.g. `"reexports.Rd"`) to leave out.
#' @param skip_internal whether to leave out pages with `\keyword{internal}`.
#' @returns a named list, with one element per Rd file of a function (named
#'   after the file), holding its sections (see [rd_sections()]).
#' @keywords internal
function_rd_sections <- function(pkg = ".", exclude = character(0),
                                 skip_internal = FALSE) {
  man_dir <- file.path(pkg, "man")
  if (!dir.exists(man_dir)) {
    stop("there is no man directory in ", normalizePath(pkg), call. = FALSE)
  }
  # only the Rd files at the top of man/ (not man/figures, macros, etc.)
  rd_files <- list.files(man_dir, pattern = "\\.[Rr]d$", full.names = TRUE)
  sections <- lapply(rd_files, rd_sections)
  names(sections) <- basename(rd_files)
  not_functions <- vapply(sections, function(s) {
    !("\\usage" %in% s) ||
      any(c("\\docType{data}", "\\docType{package}", "\\docType{import}") %in% s) ||
      (skip_internal && "\\keyword{internal}" %in% s)
  }, logical(1))
  sections[!not_functions & !(names(sections) %in% exclude)]
}

#' Report the Rd files that miss a section
#'
#' @param missing the names of the Rd files that miss the section.
#' @param section the section, as written in roxygen (e.g. `"@returns"`).
#' @param n_checked the number of Rd files that were checked.
#' @param error whether to stop with an error if some files miss the section
#'   (otherwise a warning is given).
#' @returns `missing`, invisibly.
#' @keywords internal
report_missing_section <- function(missing, section, n_checked, error) {
  if (length(missing) == 0) {
    message(
      "Success: all ", n_checked, " help pages of functions have ",
      section, "."
    )
  } else {
    msg <- paste0(
      length(missing), " of ", n_checked, " help pages of functions miss ",
      section, ":\n", paste0("  - ", missing, collapse = "\n")
    )
    if (error) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(missing)
}

#' Check that every help page of a function documents its return value
#'
#' CRAN asks that the help page of every exported function says what it
#' returns, in a `\value` section (written with `@returns` in roxygen), even
#' for functions called for their side effects (where it can say, e.g., that
#' nothing is returned, or that `x` is returned invisibly). R CMD check does
#' not check this.
#'
#' Pages of datasets, the package page and re-exported functions are not
#' checked (see [function_rd_sections()]). Run `devtools::document()` first, so
#' that the Rd files are up to date.
#'
#' @param pkg the path to the root of the package. Default is the working
#'   directory.
#' @param exclude the names of Rd files (e.g. `"reexports.Rd"`) not to check.
#' @param skip_internal whether to leave out pages with `\keyword{internal}`
#'   (roxygen `@keywords internal`). CRAN asks for these sections on the
#'   pages of exported functions; internal pages usually document functions
#'   that users do not call, but some document exported functions. Default
#'   is `FALSE` (all pages are checked).
#' @param error whether to stop with an error if some pages have no `\value`
#'   section. If `FALSE`, a warning is given instead. Default is `TRUE`.
#' @returns the names of the Rd files without a `\value` section (a character
#'   vector, empty if there are none), invisibly. A message reports the
#'   success of the check.
#' @examples
#' \dontrun{
#' check_returns_in_documentation("~/git/mypackage")
#' }
check_returns_in_documentation <- function(pkg = ".", exclude = character(0),
                                           skip_internal = FALSE,
                                           error = TRUE) {
  sections <- function_rd_sections(pkg, exclude, skip_internal)
  missing <- names(sections)[!vapply(
    sections, function(s) "\\value" %in% s, logical(1)
  )]
  report_missing_section(missing, "@returns (\\value)", length(sections), error)
}

#' Check that every help page of a function has examples
#'
#' CRAN asks that the help page of every exported function has examples, in an
#' `\examples` section (written with `@examples` in roxygen). Examples that
#' can not run on CRAN (e.g. because they need external software or take a
#' long time) should be wrapped in `\donttest{}`, or use `@examplesIf`, rather
#' than being left out. R CMD check does not check this.
#'
#' Pages of datasets, the package page and re-exported functions are not
#' checked (see [function_rd_sections()]). Run `devtools::document()` first, so
#' that the Rd files are up to date.
#'
#' @inheritParams check_returns_in_documentation
#' @param error whether to stop with an error if some pages have no
#'   `\examples` section. If `FALSE`, a warning is given instead. Default is
#'   `TRUE`.
#' @returns the names of the Rd files without an `\examples` section (a
#'   character vector, empty if there are none), invisibly. A message reports
#'   the success of the check.
#' @examples
#' \dontrun{
#' check_examples_in_documentation("~/git/mypackage")
#' }
check_examples_in_documentation <- function(pkg = ".", exclude = character(0),
                                            skip_internal = FALSE,
                                            error = TRUE) {
  sections <- function_rd_sections(pkg, exclude, skip_internal)
  missing <- names(sections)[!vapply(
    sections, function(s) "\\examples" %in% s, logical(1)
  )]
  report_missing_section(
    missing, "@examples (\\examples)", length(sections), error
  )
}
