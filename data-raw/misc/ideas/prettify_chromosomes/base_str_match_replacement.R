#' Extract all matches and capture groups
#'
#' This function is a drop-in replacement for `stringr::str_match_all()`, which
#' extracts all matches and capture groups from a character vector
#'
#' @param string A character input vector
#' @param pattern A regular expression
#' @param ignore_case Logical, should case be ignored?
#' @return A list of matrices, each matrix containing the full matches and
#' capture groups for each element of the input vector
#' @keywords internal
#' @noRd
base_str_match_all <- function(string, pattern, ignore_case = FALSE) {
  # ensure pattern uses perl regex for full regex features
  matches <- gregexpr(pattern, string, perl = TRUE, ignore.case = ignore_case)

  # extract full matches and capture groups
  lapply(seq_along(matches), function(i) {
    m <- matches[[i]]
    if (m[1] == -1) {
      return(matrix(NA_character_, 0, 1))
    } # no matches

    # extract substring positions
    starts <- as.vector(m)
    match_lens <- attr(m, "match.length")

    # get text of full matches
    full_matches <- substring(string[i], starts, starts + match_lens - 1)

    # use regexec on each full match to get capture groups
    group_list <- lapply(full_matches, function(x) {
      m2 <- regexec(pattern, x, perl = TRUE, ignore.case = ignore_case)
      regmatches(x, m2)[[1]]
    })

    # turn into matrix
    do.call(rbind, group_list)
  })
}
