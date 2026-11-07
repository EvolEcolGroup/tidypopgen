# Check the URLs of an R package as CRAN does (with urlchecker), and sort the
# problems into those that need fixing, those that need a look, and the
# usual false positives, which can be ignored.
#
# Usage, from any working directory:
#   source("<skill directory>/scripts/url_checks.R")
#   urls <- check_urls("path/to/package")
#   urls[urls$action == "fix", ]

#' The host of a URL
#'
#' @param url a character vector of URLs.
#' @returns the host of each URL, in lower case (e.g. `"github.com"`), or
#'   `""` if it can not be found.
#' @keywords internal
url_host <- function(url) {
  host <- sub("^[A-Za-z][A-Za-z0-9+.-]*://([^/?#:]+).*$", "\\1", url)
  host[host == url] <- ""
  tolower(host)
}

#' Whether hosts belong to domains
#'
#' @param host a character vector of hosts.
#' @param domains a character vector of domains (e.g. `"github.com"`); a host
#'   belongs to a domain if it is the domain or one of its subdomains (e.g.
#'   `"raw.githubusercontent.com"` belongs to `"githubusercontent.com"`).
#' @returns a logical vector, one element per host.
#' @keywords internal
host_in_domains <- function(host, domains) {
  vapply(host, function(h) {
    any(h == domains | endsWith(h, paste0(".", domains)))
  }, logical(1), USE.NAMES = FALSE)
}

#' Check the URLs of a package, and sort out the false positives
#'
#' Runs [urlchecker::url_check()], which checks the URLs of the package (in
#' DESCRIPTION, the Rd files, the vignettes, README and NEWS) as the
#' incoming checks of CRAN do, and classifies each problem:
#'
#' - `"fix"`: the URL must be changed before submitting. These are URLs that
#'   are not found (404, 410), CRAN URLs not in canonical form, and
#'   permanent redirects (301, 308) of URLs other than DOIs, which
#'   `urlchecker::url_update()` can fix.
#' - `"ignore"`: a usual false positive, which CRAN also accepts:
#'   - time outs and connection errors on well-known domains (see
#'     `reliable_domains`), which are temporary;
#'   - redirects of DOIs (`doi.org`), which always redirect to the publisher
#'     (and should not be replaced by the publisher's URL).
#' - `"review"`: anything else, to be checked by hand, e.g. 403 (Forbidden),
#'   429 (Too many requests) or 503, which publishers and some websites return
#'   to automated requests while the page works in a browser; time outs on
#'   other domains; temporary redirects (302, 307).
#'
#' @param pkg the path to the root of the package. Default is the working
#'   directory.
#' @param reliable_domains the domains on which time outs and connection
#'   errors are ignored. Default is CRAN (`r-project.org`), GitHub
#'   (`github.com`, `github.io`, `githubusercontent.com`) and Bioconductor.
#' @param doi_domains the domains of DOI resolvers, whose redirects are
#'   ignored. Default is `doi.org` (which includes `dx.doi.org`).
#' @param progress whether to show the progress of the checks.
#' @returns a data frame with one row per URL with a problem, and columns
#'   `url`, `from` (the files where the URL is, separated by commas),
#'   `status` (the HTTP status, or `"Error"`), `message`, `new` (the URL it
#'   redirects to, if any) and `action` (`"fix"`, `"review"` or `"ignore"`),
#'   ordered by action, invisibly. A message reports the number of URLs in
#'   each class, or the success of the check if there is nothing to fix or
#'   review.
#' @examples
#' \dontrun{
#' urls <- check_urls("~/git/mypackage")
#' urls[urls$action != "ignore", ]
#' }
check_urls <- function(pkg = ".",
                       reliable_domains = c(
                         "r-project.org", "github.com", "github.io",
                         "githubusercontent.com", "bioconductor.org"
                       ),
                       doi_domains = "doi.org",
                       progress = FALSE) {
  if (!requireNamespace("urlchecker", quietly = TRUE)) {
    stop("the urlchecker package is needed: install.packages(\"urlchecker\")",
      call. = FALSE
    )
  }
  res <- urlchecker::url_check(pkg, progress = progress, fail = FALSE)
  n <- NROW(res)
  urls <- data.frame(
    url = as.character(res$URL),
    from = vapply(res$From, paste, character(1), collapse = ", "),
    status = as.character(res$Status),
    message = as.character(res$Message),
    new = as.character(res$New),
    action = rep("review", n),
    stringsAsFactors = FALSE
  )
  if (n > 0) {
    # CRAN URLs not in canonical form are reported in their own column
    cran <- as.character(res$CRAN)
    non_canonical <- !is.na(cran) & nzchar(cran)
    urls$url[non_canonical] <- cran[non_canonical]
    urls$status[non_canonical] <- "Error"
    urls$message[non_canonical] <- "CRAN URL not in canonical form"

    host <- url_host(urls$url)
    status <- urls$status
    is_redirect <- status %in% c("301", "302", "303", "307", "308") |
      (!is.na(urls$new) & nzchar(urls$new))
    is_error <- !grepl("^[0-9]+$", status)
    is_doi <- host_in_domains(host, doi_domains)

    # ignored: DOI redirects, and errors (time outs, failed connections) on
    # reliable domains
    urls$action[is_doi & is_redirect] <- "ignore"
    urls$action[is_error & host_in_domains(host, reliable_domains) &
      !non_canonical] <- "ignore"
    # to fix: missing pages, non canonical CRAN URLs and permanent redirects
    # (other than DOIs)
    urls$action[status %in% c("404", "410")] <- "fix"
    urls$action[non_canonical] <- "fix"
    urls$action[status %in% c("301", "308") & !is_doi] <- "fix"
    urls <- urls[order(match(urls$action, c("fix", "review", "ignore"))), ]
    rownames(urls) <- NULL
  }

  counts <- table(factor(urls$action, levels = c("fix", "review", "ignore")))
  if (counts[["fix"]] + counts[["review"]] == 0) {
    message(
      "Success: no URL needs fixing",
      if (counts[["ignore"]] > 0) {
        paste0(" (", counts[["ignore"]], " false positives ignored)")
      },
      "."
    )
  } else {
    message(
      counts[["fix"]], " URLs to fix, ", counts[["review"]],
      " to review, and ", counts[["ignore"]], " false positives ignored."
    )
  }
  invisible(urls)
}
