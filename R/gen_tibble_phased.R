#' Constructor for a `gen_tibble` of phased genotypes
#'
#' A phased `gen_tibble` is a [`gen_tibble`] of diploid individuals that also
#' knows which allele of each genotype is on the first and which on the second
#' haplotype (e.g. `0|1` versus `1|0`), as needed by haplotype based statistics
#' (such as those of the package `tidygensel`). Phased genotypes are stored in
#' the same file backed matrix as unphased ones, with the encoding described in
#' [phased_encoding], so that all other `tidypopgen` functions (e.g.
#' [loci_alt_freq()], [loci_missingness()], [pop_fst()] or
#' [gt_pca_partialSVD()]) work on phased `gen_tibble`s as usual, treating their
#' genotypes as ordinary dosages of the alternate allele. [is_phased()] tests
#' whether a `gen_tibble` is phased, and [show_haplotypes()] returns its
#' haplotypes.
#'
#' [gen_tibble()] ignores phasing (a phased *VCF* read with [gen_tibble()]
#' gives an ordinary, unphased `gen_tibble`); use `gen_tibble_phased()` to keep
#' the phase.
#'
#' - *VCF* files: as for [gen_tibble()], only biallelic SNPs are read. All
#' individuals have to be diploid, and all genotypes phased (with `|`, e.g.
#' `0|1`): any unphased genotype (with `/`, including homozygotes such as
#' `0/0` and missing genotypes such as `./.`) is an error, naming the
#' individual and the locus. A genotype with a missing allele (e.g. `.|.`,
#' `0|.`, or a single `.`) is missing.
#'
#' - a matrix of haplotypes (mostly useful for small datasets and tests), with
#' two rows per individual (the first and second haplotypes of individual `i`
#' are rows `2 * i - 1` and `2 * i`) and one column per locus, containing 0
#' (reference allele), 1 (alternate allele) or `NA` (missing). If either
#' haplotype of an individual is missing at a locus, the whole genotype is
#' missing.
#'
#' Functions that write genotypes keep the phase where they can:
#' [gt_as_vcf()] writes phased genotypes (`0|1`), subsetting, [gt_save()],
#' [gt_load()], [gt_order_loci()] and [gt_update_backingfile()] keep phased
#' data phased, and [rbind.gen_tbl()] keeps the phase when both `gen_tibble`s
#' are phased (it then never swaps or flips alleles: it stops if the alleles
#' of a locus differ between the two datasets), and discards it (with a
#' message) when only one of them is phased. Imputation
#' ([gt_impute_simple()], [gt_impute_xgboost()]) is not available for phased
#' `gen_tibble`s; read the data with [gen_tibble()] to impute them.
#'
#' @param x can be:
#' - a string giving the path to a *VCF* file (`.vcf` or `.vcf.gz`) of phased
#'   genotypes;
#' - a matrix of haplotypes (see Details).
#' @param indiv_meta a list, data.frame or tibble with information on each
#'   individual, with a compulsory column 'id' (as for [gen_tibble()]). Only
#'   used if `x` is a matrix of haplotypes.
#' @param loci a data.frame or tibble with information on each locus, with
#'   compulsory columns 'name', 'chromosome', 'position', 'genetic_dist',
#'   'allele_ref' and 'allele_alt' (as for [gen_tibble()]). Only used if `x` is
#'   a matrix of haplotypes.
#' @param ... if `x` is the path to a *VCF* file and `parser = "vcfR"`,
#'   additional arguments passed to [vcfR::read.vcfR()]. Otherwise, unused.
#' @param parser the name of the parser used for *VCF* files, either "cpp" to
#'   use a fast C++ parser (the default), or "vcfR" to use the R package
#'   `vcfR` (see [gen_tibble()]).
#' @param chunk_size the number of loci processed at a time (currently only
#'   used with parser "vcfR").
#' @param valid_alleles a vector of valid allele values; it defaults to 'A',
#'   'T', 'C' and 'G'.
#' @param missing_alleles a vector of values in the loci table that indicate a
#'   missing value for the allele value. It defaults to '0' and '.'.
#' @param backingfile the path, including the file name without extension, for
#'   the backing files used to store the data (see [gen_tibble()]). If `x` is a
#'   *VCF* and `backingfile` is NULL, the backing files are saved next to the
#'   *VCF* file; if `x` is a matrix and `backingfile` is NULL, temporary files
#'   are used (which R deletes at the end of the session).
#' @param allow_duplicates logical. If TRUE, duplicated loci (same chromosome
#'   and position, or same name) are allowed. Default is FALSE.
#' @param names_as_int logical. If TRUE, locus names are stored as integers.
#'   Only used when `x` is a *VCF* file with parser "cpp". Default is FALSE.
#' @param quiet if TRUE, no information on the files used to store the data is
#'   printed.
#' @returns a phased [`gen_tibble`] (an object of class `gen_tbl`, whose
#'   genotypes use the code [code_phased]).
#' @seealso [phased_encoding], [is_phased()], [show_haplotypes()]
#' @rdname gen_tibble_phased
#' @export
#' @examples
#' \dontshow{
#' data.table::setDTthreads(2)
#' }
#' # Create a phased gen_tibble from a phased vcf file
#' vcf_path <- system.file("extdata", "phased_example.vcf",
#'   package = "tidypopgen"
#' )
#' x <- gen_tibble_phased(vcf_path,
#'   backingfile = tempfile("phased_"),
#'   quiet = TRUE
#' )
#' is_phased(x)
#' show_haplotypes(x)
#' # all other functions use the genotypes as usual
#' loci_alt_freq(x)
#'
#' # Create a phased gen_tibble from a matrix of haplotypes
#' haps <- matrix(c(
#'   0, 1, 1,
#'   1, 1, 0,
#'   0, 0, 1,
#'   0, 1, 1
#' ), nrow = 4, byrow = TRUE)
#' loci <- data.frame(
#'   name = paste0("rs", 1:3), chromosome = 1, position = c(10, 20, 30),
#'   genetic_dist = 0, allele_ref = "A", allele_alt = "T"
#' )
#' indiv <- data.frame(id = c("a", "b"), population = "pop1")
#' x <- gen_tibble_phased(haps, indiv_meta = indiv, loci = loci, quiet = TRUE)
#' show_haplotypes(x)
#' show_genotypes(x)
gen_tibble_phased <- function(
  x,
  ...,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  backingfile = NULL,
  allow_duplicates = FALSE,
  quiet = FALSE
) {
  UseMethod("gen_tibble_phased", x)
}

#' @export
gen_tibble_phased.default <- function(x, ...) {
  stop(
    "'x' should be the path to a vcf file, or a numeric matrix of ",
    "haplotypes",
    call. = FALSE
  )
}

#' @export
#' @rdname gen_tibble_phased
gen_tibble_phased.character <- function(
  x,
  ...,
  parser = c("cpp", "vcfR"),
  chunk_size = NULL,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  backingfile = NULL,
  allow_duplicates = FALSE,
  names_as_int = FALSE,
  quiet = FALSE
) {
  parser <- match.arg(parser)
  if (length(x) != 1 || is.na(x) || nchar(x) == 0) {
    stop("x should be the path to a single vcf file", call. = FALSE)
  }
  # check that valid alleles does not contain zero
  if ("0" %in% valid_alleles) {
    stop(paste(
      "'0' cannot be a valid allele",
      "(it is the default missing allele value!)"
    ), call. = FALSE)
  }
  is_vcf <- grepl("\\.vcf(\\.gz)?$", x, ignore.case = TRUE)
  if (!is_vcf) {
    stop(
      "phased genotypes can only be read from vcf files (.vcf or .vcf.gz): ",
      x,
      call. = FALSE
    )
  }
  # as in gen_tibble(), do not overwrite existing backing files
  if (is.null(backingfile)) {
    backingfile <- change_duplicated_file_name(x)
  } else {
    backingfile <- change_duplicated_file_name(backingfile)
  }
  x <- normalizePath(x, mustWork = FALSE)

  x_gt <- gen_tibble_vcf(
    x = x,
    ...,
    parser = parser,
    chunk_size = chunk_size,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles,
    backingfile = backingfile,
    allow_duplicates = allow_duplicates,
    names_as_int = names_as_int,
    quiet = quiet,
    phased = TRUE
  )

  # check alleles (as in gen_tibble())
  loci <- show_loci(x_gt)
  doubles <- which(is.na(loci$allele_ref) & is.na(loci$allele_alt))
  if (length(doubles) > 0) {
    check_missing <- x_gt %>%
      select_loci(all_of(doubles)) %>%
      loci_missingness()
    if (any(check_missing < 1, na.rm = TRUE)) {
      # clean up files if we are stopping (dropping the gen_tibble first, so
      # that its backing file can be unmapped, see remove_fbm_files())
      files <- gt_get_file_names(x_gt)
      rm(x_gt)
      remove_fbm_files(files)
      stop(paste(
        "Some loci are missing both reference and alternate alleles.",
        "Genotypes are not missing at these loci.",
        "Please check the loci and genotype data."
      ))
    } else {
      warning(
        "Your data contain loci with no genotypes or allele ",
        "information. Use loci_missingness() to identify them ",
        "and select_loci() to remove them."
      )
    }
  }

  gt_save(x_gt, quiet = quiet)
  x_gt
}

#' @export
#' @rdname gen_tibble_phased
gen_tibble_phased.matrix <- function(
  x,
  indiv_meta,
  loci,
  ...,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  backingfile = NULL,
  allow_duplicates = FALSE,
  quiet = FALSE
) {
  if ("ploidy" %in% names(list(...))) {
    stop("phased gen_tibbles are diploid: 'ploidy' can not be set",
      call. = FALSE
    )
  }
  rlang::check_dots_empty()
  if (!is.numeric(x) && !is.logical(x)) {
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

#' @export
gen_tibble_phased.data.frame <- function(x, ...) {
  gen_tibble_phased(as.matrix(x), ...)
}

#' Encoding of phased genotypes in a `gen_tibble`
#'
#' Phased diploid genotypes (see [gen_tibble_phased()]) are stored in the same
#' file backed matrix (a [bigstatsr::FBM.code256]) as unphased ones, using one
#' byte per genotype:
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
#' @seealso [gen_tibble_phased()], [is_phased()], [show_haplotypes()]
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
#' Phased genotypes (see [gen_tibble_phased()]) are recognised by the
#' `code256` of their file backed matrix, which is identical to [code_phased].
#'
#' @param .x a [`gen_tibble`], or its `genotypes` column.
#' @returns `TRUE` if the genotypes are phased, `FALSE` otherwise.
#' @seealso [gen_tibble_phased()], [phased_encoding]
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' is_phased(example_gt)
#'
#' vcf_path <- system.file("extdata", "phased_example.vcf",
#'   package = "tidypopgen"
#' )
#' x <- gen_tibble_phased(vcf_path,
#'   backingfile = tempfile("phased_"),
#'   quiet = TRUE
#' )
#' is_phased(x)
is_phased <- function(.x) {
  if (inherits(.x, "data.frame")) {
    stopifnot_gen_tibble(.x)
  } else if (!inherits(.x, "vctrs_bigSNP")) {
    stop(".x should be a gen_tibble, or its genotypes column", call. = FALSE)
  }
  identical(.gt_get_fbm(.x)$code256, code_phased)
}

# stop if the genotypes are not phased
stopifnot_phased <- function(.x) {
  if (!is_phased(.x)) {
    stop(
      "the genotypes are not phased (see gen_tibble_phased())",
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
#' Extract the haplotypes (as a matrix) from a phased `gen_tibble` (see
#' [gen_tibble_phased()]).
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
#' x <- gen_tibble_phased(vcf_path,
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
