# Functions for phased gen_tibbles (created with gen_tibble(phased = TRUE)):
# the encoding of phased genotypes, is_phased() and show_haplotypes()

# create a phased gen_tibble from a matrix of haplotypes x, with two rows per
# individual (rows 2 * i - 1 and 2 * i for its first and second haplotypes)
# and one column per locus, containing 0, 1 or NA (the genotype is missing if
# either haplotype is missing); used by gen_tibble(x, phased = TRUE) on a
# matrix, with the other arguments of gen_tibble()
gen_tibble_haplotypes <- function(x, indiv_meta, loci, ploidy, valid_alleles,
                                  missing_alleles, backingfile,
                                  allow_duplicates, quiet) {
  if (any(ploidy != 2)) {
    stop("phased gen_tibbles are diploid: 'ploidy' should be 2",
      call. = FALSE
    )
  }
  if (inherits(x, "data.frame")) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || (!is.numeric(x) && !is.logical(x))) {
    stop("'x' should be a numeric matrix of haplotypes", call. = FALSE)
  }
  if (any(!x %in% c(0, 1, NA))) {
    stop("'x' should only contain 0, 1 or NA", call. = FALSE)
  }
  if (nrow(x) %% 2 != 0) {
    stop("'x' should have two rows (haplotypes) per individual",
      call. = FALSE
    )
  }
  n_ind <- nrow(x) / 2
  hap1 <- x[seq(1, by = 2, length.out = n_ind), , drop = FALSE]
  hap2 <- x[seq(2, by = 2, length.out = n_ind), , drop = FALSE]
  storage.mode(hap1) <- "integer"
  storage.mode(hap2) <- "integer"
  # create an ordinary (unphased) gen_tibble of the dosages, which validates
  # the metadata, checks the dimensions and creates the backing files
  x_gt <- gen_tibble(hap1 + hap2,
    indiv_meta = indiv_meta,
    loci = loci,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles,
    backingfile = backingfile,
    allow_duplicates = allow_duplicates,
    quiet = TRUE
  )
  # then overwrite the genotypes with the phased bytes, and set the code
  bytes <- phased_genotype_byte[1 + 2 * hap1 + hap2]
  bytes[is.na(bytes)] <- phased_missing_byte
  fbm <- .gt_get_fbm(x_gt)
  fbm[] <- as.raw(bytes)
  fbm$code256 <- code_phased
  fbm$save()
  attr(x_gt$genotypes, "fbm_md5sum") <- tools::md5sum(fbm$rds)
  gt_save(x_gt, quiet = quiet)
  x_gt
}

#' Encoding of phased genotypes in a `gen_tibble`
#'
#' Phased diploid genotypes (read with `phased = TRUE` in [gen_tibble()]) are
#' stored in the same file backed matrix (a [bigstatsr::FBM.code256]) as
#' unphased ones, using one byte per genotype:
#'
#' | byte | genotype | dosage (`code256`) |
#' |------|----------|--------------------|
#' | 0    | 0\|0     | 0                  |
#' | 1    | 0\|1     | 1                  |
#' | 2    | 1\|1     | 2                  |
#' | 3    | missing  | `NA`               |
#' | 4    | 1\|0     | 1                  |
#'
#' Bytes 0 to 3 have the same meaning as in unphased diploid data (where byte 3
#' is the missing value), so that code that reads or writes the raw bytes keeps
#' working; only the second orientation of a heterozygote (`1|0`) uses a new
#' byte. The `code256` of the matrix, `code_phased`, translates each byte into
#' a dosage, so that all functions that decode genotypes through `code256` see
#' phased data as ordinary unphased dosages. `code_phased` is also how phased
#' data are recognised ([is_phased()]): it differs from the codes used for
#' unphased ([bigsnpr::CODE_012]), imputed ([bigsnpr::CODE_IMPUTE_PRED]) and
#' dosage ([bigsnpr::CODE_DOSAGE]) data, as it is the only one that decodes
#' byte 4 as 1.
#'
#' `phased_bytes` gives, for each byte, the alleles carried by the first and
#' second haplotype of an individual (`NA` for missing, and for bytes that are
#' not used).
#'
#' @format `code_phased` is a numeric vector of length 256. `phased_bytes` is
#'   an integer matrix with 256 rows (bytes 0 to 255) and two columns (`hap1`
#'   and `hap2`).
#' @seealso [gen_tibble()], [is_phased()], [show_haplotypes()]
#' @name phased_encoding
#' @examples
#' code_phased[1:6]
#' phased_bytes[1:6, ]
NULL

#' @rdname phased_encoding
#' @export
code_phased <- c(0, 1, 2, NA, 1, rep(NA_real_, 251))

#' @rdname phased_encoding
#' @export
phased_bytes <- local({
  m <- matrix(NA_integer_,
    nrow = 256, ncol = 2,
    dimnames = list(NULL, c("hap1", "hap2"))
  )
  m[1, ] <- c(0L, 0L) # byte 0: 0|0
  m[2, ] <- c(0L, 1L) # byte 1: 0|1
  m[3, ] <- c(1L, 1L) # byte 2: 1|1
  # byte 3: missing (NA, NA)
  m[5, ] <- c(1L, 0L) # byte 4: 1|0
  m
})

# the byte that encodes each phased genotype, indexed by 1 + 2 * hap1 + hap2
# (i.e. 0|0, 0|1, 1|0, 1|1), and the byte of missing genotypes (the same as
# in unphased diploid data, max_ploidy + 1)
phased_genotype_byte <- c(0L, 1L, 4L, 2L)
phased_missing_byte <- 3L

#' Test whether a `gen_tibble` contains phased genotypes
#'
#' Phased genotypes (see `phased` in [gen_tibble()]) are recognised by the
#' `code256` of their file backed matrix, which is identical to [code_phased].
#'
#' @param .x a [`gen_tibble`], or its `genotypes` column.
#' @returns `TRUE` if the genotypes are phased, `FALSE` otherwise.
#' @seealso [gen_tibble()], [phased_encoding]
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' is_phased(example_gt)
#'
#' vcf_path <- system.file("extdata", "phased_example.vcf",
#'   package = "tidypopgen"
#' )
#' x <- gen_tibble(vcf_path,
#'   phased = TRUE,
#'   backingfile = tempfile("phased_"),
#'   quiet = TRUE
#' )
#' is_phased(x)
is_phased <- function(.x) {
  if (inherits(.x, "data.frame")) {
    # use the genotypes column, so that this also works on tibbles that have
    # lost the gen_tbl class (e.g. from as_tibble() or pick())
    stopifnot_gen_tibble(.x)
    .x <- .x$genotypes
  } else if (!inherits(.x, "vctrs_bigSNP")) {
    stop(".x should be a gen_tibble, or its genotypes column", call. = FALSE)
  }
  identical(.gt_get_fbm(.x)$code256, code_phased)
}

# stop if the genotypes are not phased
stopifnot_phased <- function(.x) {
  if (!is_phased(.x)) {
    stop(
      "the genotypes are not phased (see the argument phased of gen_tibble())",
      call. = FALSE
    )
  }
  invisible(.x)
}

# stop if not all individuals are diploid (phased genotypes are only stored
# for diploid individuals)
stopifnot_phased_ploidy <- function(ploidy, sample_names) {
  not_diploid <- which(ploidy != 2)
  if (length(not_diploid) > 0) {
    stop(
      "phased genotypes can only be read for diploid individuals, but ",
      length(not_diploid), " individual(s) are not diploid (e.g. ",
      sample_names[not_diploid[1]], ", with ploidy ",
      ploidy[not_diploid[1]], ")",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Show the haplotypes of a phased `gen_tibble`
#'
#' Extract the haplotypes (as a matrix) from a phased `gen_tibble` (read with
#' `phased = TRUE` in [gen_tibble()]).
#'
#' @param .x a phased [`gen_tibble`], or its `genotypes` column.
#' @param indiv_indices the indices of the individuals (rows of `.x`) to show.
#' @param loci_indices the indices of the loci to show.
#' @returns an integer matrix with two rows per individual (its first and
#'   second haplotypes, named `<id>_1` and `<id>_2`) and one column per locus
#'   (named after the locus), with 0 for the reference allele, 1 for the
#'   alternate allele and `NA` for missing genotypes.
#' @seealso [show_genotypes()] for the dosages of the alternate allele.
#' @export
#' @examples
#' vcf_path <- system.file("extdata", "phased_example.vcf",
#'   package = "tidypopgen"
#' )
#' x <- gen_tibble(vcf_path,
#'   phased = TRUE,
#'   backingfile = tempfile("phased_"),
#'   quiet = TRUE
#' )
#' show_haplotypes(x)
#' show_haplotypes(x, indiv_indices = 1:2, loci_indices = 2:4)
show_haplotypes <- function(.x, indiv_indices = NULL, loci_indices = NULL) {
  stopifnot_phased(.x)
  if (inherits(.x, "data.frame")) {
    genotypes <- .x$genotypes
    ids <- as.character(.x$id)
  } else {
    genotypes <- .x
    ids <- names(genotypes)
    if (is.null(ids)) {
      ids <- as.character(seq_along(genotypes))
    }
  }
  fbm <- .gt_get_fbm(genotypes)
  rows <- vctrs::vec_data(genotypes)
  loci <- show_loci(genotypes)
  cols <- loci$big_index
  loci_names <- loci$name
  if (!is.null(indiv_indices)) {
    rows <- rows[indiv_indices]
    ids <- ids[indiv_indices]
  }
  if (!is.null(loci_indices)) {
    cols <- cols[loci_indices]
    loci_names <- loci_names[loci_indices]
  }
  # read the raw bytes (a code256 of 0:255 returns the bytes themselves)
  raw_fbm <- fbm$copy(code = as.numeric(0:255))
  bytes <- raw_fbm[rows, cols, drop = FALSE]
  haps <- matrix(NA_integer_, nrow = 2 * length(rows), ncol = length(cols))
  haps[seq(1, by = 2, length.out = length(rows)), ] <-
    phased_bytes[bytes + 1, 1]
  haps[seq(2, by = 2, length.out = length(rows)), ] <-
    phased_bytes[bytes + 1, 2]
  rownames(haps) <- paste0(
    rep(ids, each = 2),
    rep(c("_1", "_2"), times = length(ids))
  )
  colnames(haps) <- as.character(loci_names)
  haps
}
