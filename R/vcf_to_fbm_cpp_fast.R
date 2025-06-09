#' Convert vcf to FBM.
#'
#' Convert a vcf file to a Filebacked Big Matrix (FBM) object. This should work
#' even for large vcf files that would not fit in memory.
#'
#' @param vcf_path the path to the vcf
#' @param backingfile the name of the file to use as the backing file for the
#'   FBM. If NULL, the vcf path will be used.
#' @return path to the resulting rds file as class bigSNP.
#' @keywords internal
#' @noRd
vcf_to_fbm_cpp <- function(
    vcf_path,
    backingfile = NULL,
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


  # individual metadata table
  fam <- tibble(
    family.ID = vcf_meta$sample_names,
    sample.ID = vcf_meta$sample_names,
    paternal.ID = 0,
    maternal.ID = 0,
    sex = 0,
    affection = -9,
    ploidy = ploidy
  )


  # loci metadata table
  loci <- vcf_meta$loci_tbl

  # add an empty genetic.pos column
  loci <- loci %>% mutate(genetic.dist = 0, .after = "physical.pos")
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

  bigsnp_obj <- structure(
    list(
      genotypes = file_backed_matrix,
      fam = fam,
      map = loci
    ),
    class = "bigSNP"
  )

  bigsnp_obj <- bigsnpr::snp_save(bigsnp_obj)
  # and return the path to the rds
  bigsnp_obj$genotypes$rds
}
