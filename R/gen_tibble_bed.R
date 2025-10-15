gen_tibble_bed_rds <- function(
    x,
    ...,
    ploidy = 2,
    valid_alleles = c("A", "T", "C", "G"),
    missing_alleles = c("0", "."),
    allow_duplicates = FALSE,
    backingfile = NULL,
    quiet = FALSE) {
  # if file ends in bed
  if (grepl("\\.bed$", x)) {
    # check that the bim and fam files exist
    if (!file.exists(bigsnpr::sub_bed(x, ".bim"))) {
      stop("bim file ", bigsnpr::sub_bed(x, ".bim"), " does not exist")
    }
    if (!file.exists(bigsnpr::sub_bed(x, ".fam"))) {
      stop("fam file ", bigsnpr::sub_bed(x, ".fam"), " does not exist")
    }
    # read in the loci and indiv_meta from the bim and fam files
    loci <- loci_from_bim(bigsnpr::sub_bed(x, ".bim"),
      valid_alleles = valid_alleles,
      missing_alleles = missing_alleles,
      allow_duplicates = allow_duplicates
    )
    indiv_meta <- indiv_meta_from_fam(bigsnpr::sub_bed(x, ".fam"))

    if (is.null(backingfile)) {
      backingfile <- bigsnpr::sub_bed(x)
    }
    fbm_path <- fbm_read_bed(x,
      n_indiv = nrow(indiv_meta),
      n_snp = nrow(loci),
      backingfile = backingfile
    )

    fbm_obj <- readRDS(fbm_path)
  } else {
    # read the bigsnp object from the rds file
    bigsnp_obj <- bigsnpr::snp_attach(x)
    fbm_obj <- bigsnp_obj$genotypes
    # note that position and genetic dist are already flipped in the bigsnpr obj
    loci <- loci_from_bim(bigsnp_obj$map,
      valid_alleles = valid_alleles,
      missing_alleles = missing_alleles,
      allow_duplicates = allow_duplicates
    )
    indiv_meta <- indiv_meta_from_fam(bigsnp_obj$fam)
    # create a copy of the bigsnp object
    # new file name for the bigsnp ends in "_bigsnp.rds"
    file.copy(x, sub("\\.bk$", "_bigsnp.rds", fbm_obj$backingfile))
    fbm_path <- sub("\\.bk$", ".rds", fbm_obj$backingfile)
    # replace the bigsnp rds with the fbm rds
    saveRDS(fbm_obj, fbm_path)
    # remove the bisnp object from memory
    rm(bigsnp_obj)
  }
  # add ploidy check and fbm_ploidy
  if (length(ploidy) > 1) {
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
  } else {
    if ((ploidy != 0) && (ploidy != -2)) {
      fbm_ploidy <- rep(ploidy, nrow(indiv_meta))
    } else {
      stop(
        "'ploidy' 0 (mixed ploidy) or -2 (haplodiploids) ",
        "require a vector of individual ploidies"
      )
    }
  }

  indiv_meta$genotypes <- new_vctrs_bigsnp(
    fbm_obj = fbm_obj,
    fbm_file = fbm_path,
    loci = loci,
    indiv_id = indiv_meta$id,
    ploidy = ploidy,
    fbm_ploidy = fbm_ploidy
  )

  new_gen_tbl <- tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )
  return(new_gen_tbl)
}

# this is a stripped down version of bigsnpr::snp_readBed
# it only creates the fbm of the genotypes
fbm_read_bed <- function(bedfile,
                         n_indiv,
                         n_snp,
                         backingfile = bigsnpr::sub_bed(bedfile)) {
  backingfile <- path.expand(backingfile)
  # check that file exists
  if (file.exists(paste0(backingfile, ".bk"))) {
    stop("the backing file ", paste0(backingfile, ".bk"), " already exists")
  }
  if (file.exists(paste0(backingfile, ".rds"))) {
    stop("the backing file ", paste0(backingfile, ".rds"), " already exists")
  }
  big_geno <- bigstatsr::FBM.code256(
    nrow = n_indiv, ncol = n_snp,
    code = bigsnpr::CODE_012, backingfile = backingfile, init = NULL,
    create_bk = TRUE
  )
  reach_eof <- readbina(path.expand(bedfile), big_geno, getCode())
  if (!reach_eof) {
    warning("EOF of bedfile has not been reached.")
  }
  rds <- bigstatsr::sub_bk(big_geno$backingfile, ".rds")
  saveRDS(big_geno, rds)
  rds
}








#' Generate a tibble of loci from a PLINK .bim file
#'
#' This function reads a PLINK .bim file and extracts the relevant information
#' to create a loci tibble compatible with a `gen_tibble` object.
#' @param bim the path to a bim file or a data.frame containing the contents of
#'   a bim file.
#' @param valid_alleles a character vector of valid alleles. Default is c("A",
#'   "T", "C", "G").
#' @param missing_alleles a character vector of alleles to be treated as
#'   missing. Default is c("0", ".").
#' @param allow_duplicates logical, if FALSE (default) an error is raised if
#'   there are duplicate locus names. If TRUE, duplicate locus names are
#'   allowed, but a warning is issued.
#' @returns A tibble with columns: `big_index`, `name`, `chromosome`,
#'   `position`, `gentic_dist`, `allele_ref`, `allele_alt`
#' @keywords internal
#' @noRd

loci_from_bim <- function(bim,
                          valid_alleles = c("A", "T", "C", "G"),
                          missing_alleles = c("0", "."),
                          allow_duplicates = FALSE) {
  if (is.character(bim)) {
    bim <- utils::read.table(bim, stringsAsFactors = FALSE)
    names(bim)[1:6] <- c(
      "chromosome", "marker.ID", "genetic.dist",
      "physical.pos", "allele1", "allele2"
    )
  }
  if (!is.data.frame(bim) || ncol(bim) < 6) {
    stop(
      "bim must be a data.frame with at least 6 columns or a path",
      "to a .bim file"
    )
  }
  loci <- tibble::tibble(
    name = as.character(bim$marker.ID),
    chromosome = as.factor(bim$chromosome),
    position = as.integer(bim$physical.pos),
    genetic_dist = as.numeric(bim$genetic.dist),
    allele_ref = as.character(bim$allele2),
    allele_alt = as.character(bim$allele1)
  )

  loci <- validate_loci(loci,
    check_alphabet = TRUE,
    check_duplicates = TRUE,
    allow_duplicates = allow_duplicates,
    harmonise_loci = TRUE,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles
  )

  return(loci)
}


#' Generate a tibble of indiv_meta from a fam file
#'
#' This function reads a PLINK .fam file and extracts the relevant information
#' to create an indiv_meta tibble compatible with a `gen_tibble` object.
#' @param fam the path to a fam file or a data.frame containing the contents of
#' a fam file.
#' @returns A tibble with columns: `id`, `population`, `paternal_ID`,
#' `maternal_ID`, `sex`, `phenotype`
#' @keywords internal
#' @noRd

indiv_meta_from_fam <- function(fam) {
  if (is.character(fam)) {
    fam <- utils::read.table(fam, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(fam) || ncol(fam) < 6) {
    stop(
      "fam must be a data.frame with at least 6 columns or a path",
      "to a .fam file"
    )
  }

  # set column names
  names(fam)[1:6] <- c(
    "family_id", "sample_id", "paternal_id", "maternal_id", "sex",
    "affection"
  )

  # start building the indiv_meta list
  indiv_meta <- list(
    id = as.character(fam$sample_id)
  )
  # if sample_id and family_id are different, add population
  if (!all(fam$family_id == fam$sample_id)) {
    indiv_meta$population <- as.character(fam$family_id)
  }

  # transfer some of the fam info to the metadata table if it is not missing
  # (0 is the default missing value)

  # if paternal_id is not 0, add paternal_ID
  if (!all(fam$paternal_id == 0)) {
    indiv_meta$paternal_ID <- as.character(fam$paternal_id)
    indiv_meta$paternal_ID[indiv_meta$paternal_ID == "0"] <-
      NA
  }
  # if maternal_id is not 0, add maternal_ID
  if (!all(fam$maternal_id == 0)) {
    indiv_meta$maternal_ID <- as.character(fam$maternal_id)
    indiv_meta$maternal_ID[indiv_meta$maternal_ID == "0"] <-
      NA
  }

  # if sex is numeric or integer
  if (inherits(fam$sex, "numeric") || inherits(fam$sex, "integer")) {
    if (!all(fam$sex == 0)) {
      indiv_meta$sex <- dplyr::case_match(
        fam$sex,
        1 ~ "male",
        2 ~ "female",
        .default = NA,
        .ptype = factor(levels = c("female", "male"))
      )
    }
  }

  if (inherits(fam$affection, "numeric") || inherits(fam$sex, "integer")) {
    if (!all(fam$affection %in% c(0, -9))) {
      indiv_meta$phenotype <- dplyr::case_match(
        fam$affection,
        1 ~ "control",
        2 ~ "case",
        -9 ~ NA,
        .default = NA,
        .ptype = factor(levels = c("control", "case"))
      )
    }
  }

  indiv_meta <- validate_indiv_meta(as.data.frame(indiv_meta))

  return(indiv_meta)
}
