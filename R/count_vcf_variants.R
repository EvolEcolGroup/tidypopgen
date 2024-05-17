#' Counts the number of variants in a vcf file
#'
#' Count the number of VCF variants by first parsing the header, and then
#' rapidly counting the number of platform-independent
#' newlines (CR, LF, and CR+LF), including a last line with neither.
#'
#' @author author Henrik Bengtsson for the original `countLines` in `R.utils`;
#' Andrea Manica for the modified version focussed on vcf
#' @param file name and path of vcf file (it can be compressed)
#' @param chunk_size The number of bytes read in each chunk.
#' @returns the number of variants (an integer)
#' @keywords internal

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

  LF <- as.raw(0x0a)
  CR <- as.raw(0x0d)
  SPC <- as.raw(32L)

  isLastCR <- isLastLF <- FALSE
  isEmpty <- TRUE
  nbrOfLines <- 0L
  while (TRUE) {
    bfr <- readBin(con = con, what = raw(), n = chunk_size)
    if (isLastCR) {
      # Don't count LF following a CR in previous chunk.
      if (bfr[1L] == LF) {
        bfr[1L] <- SPC
      }
    }

    n <- length(bfr)
    if (n == 0L) {
      break
    }

    isEmpty <- FALSE

    # Replace all CRLF:s to become LF:s
    idxsCR <- which(bfr == CR)
    nCR <- length(idxsCR)
    if (nCR > 0L) {
      idxsCRLF <- idxsCR[(bfr[idxsCR + 1L] == LF)]
      if (length(idxsCRLF) > 0L) {
        bfr <- bfr[-idxsCRLF]
        n <- length(bfr)
        idxsCRLF <- NULL # Not needed anymore
        nCR <- length(which(bfr == CR))
      }
    }

    # Count all CR:s and LF:s
    nLF <- length(which(bfr == LF))
    nbrOfLines <- nbrOfLines + (nCR + nLF)

    if (n == 0L) {
      isLastCR <- isLastLF <- FALSE
    } else {
      # If last symbol is CR it might be followed by a LF in
      # the next chunk. If so, don't count that next LF.
      bfrN <- bfr[n]
      isLastCR <- (bfrN == CR)
      isLastLF <- (bfrN == LF)
    }
  } # while()

  # Count any last line without newline too
  if (!isEmpty) {
    if (!isLastLF) nbrOfLines <- nbrOfLines + 1L
    attr(nbrOfLines, "lastLineHasNewline") <- isLastLF
  }

  # note the +1, since we read one variant whilst checking for the header
  nbrOfLines + 1
}


#' Counts the number of individuals in a vcf file
#'
#' Count the number of VCF individuals by first parsing the header, and then
#' getting the first line, then separating on tabs.
#'
#' @author author Henrik Bengtsson for the original `countLines` in `R.utils`;
#' Andrea Manica for the modified version focussed on vcf; Max Carter-Brown modified
#' the file above.
#' @param file name and path of vcf file (it can be compressed)
#' @returns the number of individuals (an integer)
#' @keywords internal

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
