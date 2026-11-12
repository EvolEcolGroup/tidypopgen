#' Convert vcf to FBM.
#'
#' Convert a vcf file to a Filebacked Big Matrix (FBM) object. This should work
#' even for large vcf files that would not fit in memory.
#'
#' @param vcf_path the path to the vcf
#' @param backingfile the name of the file to use as the backing file for the
#'   FBM. If NULL, the vcf path will be used.
#' @param valid_alleles a character vector of valid alleles. Default is c("A",
#'   "T", "C", "G").
#' @param missing_alleles a character vector of alleles to be treated as
#'   missing. Default is c("0", ".").
#' @param allow_duplicates whether to allow duplicated loci (same chromosome and
#'   position) or duplicated locus names. Default is FALSE.
#' @param names_as_int whether to convert locus names to
#'   integers (their original names are written to the locus sidecar of the
#'   backing file, see [gt_lookup_loci()]). Default is FALSE.
#' @param chromosomes_as_int whether to store the chromosomes as integer codes
#'   (their names are written to the chromosome sidecar of the backing file).
#'   Default is FALSE.
#' @param n_cores the number of cores used to decompress a bgzipped vcf and to
#'   parse the genotypes. Default is 1.
#' @param quiet whether to print messages.
#' @param phased whether to store phased diploid genotypes (see
#'   [gen_tibble()]), rather than alternate allele counts.
#' @returns an object of the class `gen_tbl`.
#' @keywords internal
#' @noRd
vcf_to_fbm_cpp <- function(
  vcf_path,
  backingfile = NULL,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  allow_duplicates = FALSE,
  names_as_int = FALSE,
  chromosomes_as_int = FALSE,
  n_cores = 1,
  quiet = FALSE,
  phased = FALSE
) {
  if (is.null(backingfile)) {
    backingfile <- vcf_path
    backingfile <- sub("\\.vcf\\.gz$", "", backingfile, ignore.case = TRUE)
    backingfile <- sub("\\.vcf$", "", backingfile, ignore.case = TRUE)
  }
  if (file_ext(backingfile) == "bk") {
    backingfile <- bigstatsr::sub_bk(backingfile, "")
  }

  if (!is.numeric(n_cores) || length(n_cores) != 1 || !is.finite(n_cores) ||
        n_cores < 1) {
    stop("'n_cores' should be a positive integer")
  }

  # the genotypes are written to the backing file while reading the VCF
  bk_file <- path.expand(paste0(backingfile, ".bk"))
  if (file.exists(bk_file)) {
    stop("the backing file ", bk_file, " already exists")
  }
  # with integer names, the original names (with the chromosomes and
  # positions) are written to the locus sidecar of the backing file, and
  # with integer chromosomes, their names to the chromosome sidecar (see
  # R/loci_sidecar.R); sidecars already next to the new backing file are left
  # from older files with the same name, and can not belong to the new data
  loci_file <- sidecar_path_from_bk(bk_file, "loci")
  chromosomes_file <- sidecar_path_from_bk(bk_file, "chromosomes")
  remove_stale_sidecars(bk_file)
  # remove the backing file (and sidecars) if we stop before creating the FBM
  # (also if the VCF was read, but the loci table could not be created)
  fbm_created <- FALSE
  on.exit(
    if (!fbm_created) unlink(c(bk_file, loci_file, chromosomes_file)),
    add = TRUE
  )
  # read the VCF once: create a loci_table, figure out no_individuals and
  # ploidy from the first marker, and write the genotypes to the backing file
  # (locus names are not read if we store them as integers, which saves a lot
  # of memory with many loci: they are written to the sidecar instead; the
  # same goes for the names of the chromosomes with integer chromosomes).
  # With n_cores > 1, the blocks of a bgzipped VCF are decompressed in
  # parallel, and the genotypes parsed in parallel. Phased genotypes are
  # stored as the bytes of phased_genotype_byte (missing genotypes are
  # max_ploidy + 1 = 3, as in unphased diploid data)
  vcf_meta <- vcf_to_bk(vcf_path,
    bk_file = bk_file, read_names = !names_as_int,
    n_threads = n_cores,
    loci_file = if (names_as_int) loci_file else "",
    chromosomes_as_int = chromosomes_as_int,
    chromosomes_file = if (chromosomes_as_int) chromosomes_file else "",
    phased_lut = if (phased) phased_genotype_byte else integer(0)
  )
  # an FBM can not have zero loci (its empty backing file can not be mapped)
  if (length(vcf_meta$loci_tbl$physical.pos) == 0) {
    stop("the VCF has no biallelic SNPs")
  }
  # figure out no_individuals and ploidy
  ploidy <- vcf_meta$ploidy
  no_individuals <- length(ploidy)

  # check that there are no missing values in ploidy vector
  if (any(is.na(ploidy))) {
    stop("'ploidy' can not contain NAs")
  }
  # check all values > 0
  if (any(ploidy <= 0)) {
    stop(
      "the vector of individual ploidies ('ploidy') must contain ",
      "positive integers"
    )
  }
  fbm_ploidy <- ploidy
  max_ploidy <- max(ploidy)
  # phased genotypes are only stored for diploid individuals
  if (phased) {
    stopifnot_phased_ploidy(ploidy, vcf_meta$sample_names)
  }

  indiv_meta <- list(
    id = vcf_meta$sample_names
  )

  # loci metadata table
  loci <- vcf_meta$loci_tbl
  # chromosomes are given as an index into their (unique) names, in order of
  # appearance: turn them into a factor (with integer chromosomes, they are
  # already their integer codes)
  if (!chromosomes_as_int) {
    loci$chromosome <- chromosome_index_to_factor(
      loci$chromosome,
      vcf_meta$chromosome_names
    )
  }
  # missing loci names were replaced in C++ by the chromosome and position
  # (as vcfR does); integer names are the index of each locus
  if (names_as_int) {
    loci$marker.ID <- seq_along(loci$physical.pos)
  }

  # (with an empty genetic distance column)
  loci <- tibble::tibble(
    name = loci$marker.ID,
    chromosome = loci$chromosome,
    position = loci$physical.pos,
    genetic_dist = 0,
    allele_ref = loci$allele2,
    allele_alt = loci$allele1
  )
  # a missing allele "." (e.g. the ALT of a monomorphic locus) is always
  # missing, also with custom missing_alleles, as in the vcfR parser (which
  # reads it as NA)
  for (allele_col in c("allele_ref", "allele_alt")) {
    missing_dot <- which(loci[[allele_col]] == ".")
    if (length(missing_dot) > 0) {
      loci[[allele_col]][missing_dot] <- NA
    }
  }

  # validate the loci (before creating the FBM)
  loci <- validate_loci(loci,
    check_alphabet = TRUE,
    harmonise_loci = TRUE,
    check_duplicates = TRUE,
    allow_duplicates = allow_duplicates,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles,
    chromosomes_as_int = chromosomes_as_int
  )
  # validate individuals
  indiv_meta <- validate_indiv_meta(as.data.frame(indiv_meta))

  if (phased) {
    code256 <- code_phased
  } else {
    code256 <- rep(NA_real_, 256)
    code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)
  }

  # create the file backed matrix, from the backing file written above
  file_backed_matrix <- bigstatsr::FBM.code256(
    nrow = no_individuals,
    ncol = nrow(loci),
    code = code256,
    backingfile = backingfile,
    create_bk = FALSE
  )
  fbm_created <- TRUE

  # construct path, and save the FBM (before its md5 sum is taken by
  # new_vctrs_bigsnp(); this takes the same time whatever the size of the
  # data, as the .rds only holds the description of the FBM)
  fbm_path <- bigstatsr::sub_bk(file_backed_matrix$backingfile, ".rds")
  saveRDS(file_backed_matrix, fbm_path)

  # create genotypes column
  indiv_meta$genotypes <- new_vctrs_bigsnp(
    fbm_obj = file_backed_matrix,
    fbm_file = fbm_path,
    loci = loci,
    indiv_id = indiv_meta$id,
    ploidy = max_ploidy,
    fbm_ploidy = fbm_ploidy
  )

  new_gen_tbl <- tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )
  return(new_gen_tbl)
}

# the chromosomes of the loci of a VCF as a factor, from their index (1-based)
# into the unique chromosome names (in order of appearance, as returned by
# vcf_to_bk()), with the levels sorted as in as.factor() (this avoids creating
# and sorting a string for each locus). A missing chromosome "." is not a
# level, so that it becomes NA, as in the vcfR parser.
chromosome_index_to_factor <- function(index, chromosome_names) {
  level_order <- order_as_levels(chromosome_names)
  level_order <- level_order[chromosome_names[level_order] != "."]
  # the code of each chromosome name in the factor is its rank among the
  # levels (NA for ".")
  chromosome_codes <- rep(NA_integer_, length(chromosome_names))
  chromosome_codes[level_order] <- seq_along(level_order)
  structure(
    chromosome_codes[index],
    levels = chromosome_names[level_order],
    class = "factor"
  )
}

# the order of a character vector of unique strings (with no NAs) as the
# levels of as.factor(), i.e. in the collation of the locale (as order()
# does): x[order_as_levels(x)] are the levels.
# Sorting millions of strings (e.g. the chromosomes of RAD loci assembled de
# novo) in the collation of the locale is slow, so we first sort them in C
# order (method = "radix", much faster); if each string then comes strictly
# after the previous one in the collation of the locale, which only needs a
# comparison per pair of neighbours (made by is.unsorted(), without the
# copies of the vector that comparing x[-1] with x[-n] needs), this is the
# order in the collation of the locale (e.g. for names that only differ in
# digits). Otherwise (e.g. with mixed case or punctuation) we sort them in the
# collation of the locale.
order_as_levels <- function(x) {
  x_order <- order(x, method = "radix")
  if (length(x) < 2 || !is.unsorted(x[x_order], strictly = TRUE)) {
    return(x_order)
  }
  order(x)
}
