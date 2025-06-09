#' Counts the number of variants in a vcf file
#'
#' Count the number of VCF variants by first parsing the header, and then
#' rapidly counting the number of platform-independent
#' newlines (carriage_return, line_feed, and carriage_return+line_feed),
#' including a last line with neither.
#'
#' @author author Henrik Bengtsson for the original `countLines` in `R.utils`;
#' Andrea Manica for the modified version focussed on vcf
#' @param file name and path of vcf file (it can be compressed)
#' @param chunk_size The number of bytes read in each chunk.
#' @returns the number of variants (an integer)
#' @keywords internal
#' @noRd

count_vcf_variants <- function(file, chunk_size = 50e6) {
  if (!is.character(file)) {
    stop("file should be a character giving the path of the vcf")
  }
  if (!file.exists(file)) {
    stop("file is not a valid path to a vcf file")
  }
  con <- gzfile(file, open = "rb")
  on.exit(close(con))

  # skip the header (all lines start with #)
  while (TRUE) {
    a_line <- readLines(con, n = 1)
    if (!substr(a_line, 1, 1) == "#") {
      break
    }
  }

  line_feed <- as.raw(0x0a)
  carriage_return <- as.raw(0x0d)
  spc <- as.raw(32L)

  is_last_carriage_return <- is_last_line_feed <- FALSE
  is_empty <- TRUE
  nbr_of_lines <- 0L
  while (TRUE) {
    bfr <- readBin(con = con, what = raw(), n = chunk_size)
    if (is_last_carriage_return) {
      # Don't count line_feed following a carriage_return in previous chunk.
      if (bfr[1L] == line_feed) {
        bfr[1L] <- spc
      }
    }

    n <- length(bfr)
    if (n == 0L) {
      break
    }

    is_empty <- FALSE

    # Replace all CRLF:s to become line_feed:s
    idxs_carriage_return <- which(bfr == carriage_return)
    n_carriage_return <- length(idxs_carriage_return)
    if (n_carriage_return > 0L) {
      idxs_cr_lf <- idxs_carriage_return[
        (bfr[idxs_carriage_return + 1L] == line_feed)
      ] # nolint
      if (length(idxs_cr_lf) > 0L) {
        bfr <- bfr[-idxs_cr_lf]
        n <- length(bfr)
        idxs_cr_lf <- NULL # Not needed anymore
        n_carriage_return <- length(which(bfr == carriage_return))
      }
    }

    # Count all carriage_return:s and line_feed:s
    n_line_feed <- length(which(bfr == line_feed))
    nbr_of_lines <- nbr_of_lines + (n_carriage_return + n_line_feed)

    if (n == 0L) {
      is_last_carriage_return <- is_last_line_feed <- FALSE
    } else {
      # If last symbol is carriage_return it might be followed by a line_feed in
      # the next chunk. If so, don't count that next line_feed.
      bfr_n <- bfr[n]
      is_last_carriage_return <- (bfr_n == carriage_return)
      is_last_line_feed <- (bfr_n == line_feed)
    }
  } # while()

  # Count any last line without newline too
  if (!is_empty) {
    if (!is_last_line_feed) nbr_of_lines <- nbr_of_lines + 1L
    attr(nbr_of_lines, "lastLineHasNewline") <- is_last_line_feed
  }

  # note the +1, since we read one variant whilst checking for the header
  nbr_of_lines + 1
}


#' Counts the number of individuals in a vcf file
#'
#' Count the number of VCF individuals by first parsing the header, and then
#' getting the first line, then separating on tabs.
#'
#' @author author Henrik Bengtsson for the original `countLines` in `R.utils`;
#'   Andrea Manica for the modified version focussed on vcf; Max Carter-Brown
#'   modified the file above.
#' @param file name and path of vcf file (it can be compressed)
#' @returns the number of individuals (an integer)
#' @keywords internal
#' @noRd
count_vcf_individuals <- function(file) {
  if (!is.character(file)) {
    stop("file should be a character giving the path of the vcf")
  }
  if (!file.exists(file)) {
    stop("file is not a valid path to a vcf file")
  }
  con <- gzfile(file, open = "rb")
  on.exit(close(con))

  # skip the header (all lines start with #)
  while (TRUE) {
    a_line <- readLines(con, n = 1)
    if (!substr(a_line, 1, 1) == "#") {
      break
    }
  }

  # read the first variant line
  a_line <- readLines(con, n = 1)
  # split the line by tab
  a_line <- strsplit(a_line, "\t")[[1]]
  # count the number of individuals
  sum(a_line[10:length(a_line)] != ".")
}
