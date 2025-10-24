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
#' @param quiet whether to print messages.
#' @returns an object of the class `gen_tbl`.
#' @keywords internal
#' @noRd
vcf_to_fbm_cpp <- function(
    vcf_path,
    backingfile = NULL,
    valid_alleles = c("A", "T", "C", "G"),
    missing_alleles = c("0", "."),
    allow_duplicates = FALSE,
    quiet = FALSE) {
  if (is.null(backingfile)) {
    backingfile <- vcf_path
    backingfile <- sub("\\.vcf.gz$", "", backingfile)
    backingfile <- sub("\\.vcf$", "", backingfile)
  }
  if (file_ext(backingfile) == "bk") {
    backingfile <- bigstatsr::sub_bk(backingfile, "")
  }

  # create a loci_table, and figure out no_individuals and ploidy from first
  # marker
  vcf_meta <- vcf_loci_table(vcf_path)
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

  code256 <- rep(NA_real_, 256)
  code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)

  # create the file backed matrix
  file_backed_matrix <- bigstatsr::FBM.code256(
    nrow = no_individuals,
    ncol = nrow(vcf_meta$loci_tbl),
    code = code256,
    backingfile = backingfile
  )

  vcf_genotypes_to_fbm(vcf_path, file_backed_matrix,
    biallelic = vcf_meta$biallelic,
    n_header_lines = vcf_meta$n_header_lines,
    missing_value = max_ploidy + 1
  )

  indiv_meta <- list(
    id = vcf_meta$sample_names
  )

  # loci metadata table
  loci <- vcf_meta$loci_tbl

  # add an empty genetic.pos column
  loci <- loci %>% mutate(genetic.dist = 0, .before = "physical.pos")
  loci$physical.pos <- as.integer(loci$physical.pos)
  # if loci names are missing, create a name with scaffold and position
  # (same behaviour as vcfR)
  missing_loci_ids <- which(loci$marker.ID == ".")
  if (length(missing_loci_ids) > 0) {
    loci$marker.ID[missing_loci_ids] <- paste(
      loci$chromosome[missing_loci_ids],
      loci$physical.pos,
      sep = "_"
    )
  }

  loci <- tibble::tibble(
    name = loci$marker.ID,
    chromosome = loci$chromosome,
    position = loci$physical.pos,
    genetic_dist = loci$genetic.dist,
    allele_ref = loci$allele2,
    allele_alt = loci$allele1
  )

  # validate the loci
  loci <- validate_loci(loci,
    check_alphabet = TRUE,
    harmonise_loci = TRUE,
    check_duplicates = TRUE,
    allow_duplicates = allow_duplicates,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles
  )
  # validate individuals
  indiv_meta <- validate_indiv_meta(as.data.frame(indiv_meta))
  # construct path
  fbm_path <- bigstatsr::sub_bk(file_backed_matrix$backingfile, ".rds")

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
