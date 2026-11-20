gen_tibble_bed_rds <- function(
  x,
  ...,
  ploidy = 2,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  allow_duplicates = FALSE,
  backingfile = NULL,
  names_as_int = FALSE,
  n_cores = 1,
  quiet = FALSE
) {
  # if file ends in bed
  if (grepl("\\.bed$", x, ignore.case = TRUE)) {
    # check that the bim and fam files exist
    bim_file <- sub_ext_ignore_case(x, "bed", ".bim")
    fam_file <- sub_ext_ignore_case(x, "bed", ".fam")
    if (!file.exists(bim_file)) {
      stop("bim file ", bim_file, " does not exist")
    }
    if (!file.exists(fam_file)) {
      stop("fam file ", fam_file, " does not exist")
    }
    # read in the loci and indiv_meta from the bim and fam files
    loci <- loci_from_bim(bim_file,
      valid_alleles = valid_alleles,
      missing_alleles = missing_alleles,
      allow_duplicates = allow_duplicates,
      names_as_int = names_as_int
    )
    indiv_meta <- indiv_meta_from_fam(fam_file)

    if (is.null(backingfile)) {
      backingfile <- sub("\\.bed$", "", x, ignore.case = TRUE)
    }
    fbm_path <- fbm_read_bed(x,
      n_indiv = nrow(indiv_meta),
      n_snp = nrow(loci),
      backingfile = backingfile,
      n_cores = n_cores
    )

    fbm_obj <- readRDS(fbm_path)
  } else {
    # read the bigsnp object from the rds file
    bigsnp_obj <- bigsnpr::snp_attach(x)
    # note that position and genetic dist are already flipped in the bigsnpr obj
    loci <- loci_from_bim(bigsnp_obj$map,
      valid_alleles = valid_alleles,
      missing_alleles = missing_alleles,
      allow_duplicates = allow_duplicates,
      names_as_int = names_as_int
    )
    indiv_meta <- indiv_meta_from_fam(bigsnp_obj$fam)
    if (is.null(backingfile)) {
      backingfile <- change_duplicated_file_name(x)
    }
    # copy the backing file of the bigsnp object, so that the fbm rds does
    # not overwrite the bigsnp rds (an FBM rds has the name of its .bk file)
    # the copy belongs to the gen_tibble, so it does not keep the permissions
    # of the original (which might be read-only), as we need to write into it
    # (e.g. to impute)
    if (!file.copy(
      bigsnp_obj$genotypes$backingfile,
      paste0(backingfile, ".bk"),
      copy.mode = FALSE
    )) {
      stop("could not copy the backing file to ", paste0(backingfile, ".bk"))
    }
    fbm_obj <- bigstatsr::FBM.code256(
      nrow = nrow(bigsnp_obj$genotypes),
      ncol = ncol(bigsnp_obj$genotypes),
      code = bigsnp_obj$genotypes$code256,
      backingfile = backingfile,
      create_bk = FALSE
    )
    fbm_path <- bigstatsr::sub_bk(fbm_obj$backingfile, ".rds")
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
        "'ploidy' 0 (mixed ploidy) or -2 (pseudohaploids) ",
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

# Read the genotypes of a PLINK .bed file into a new FBM.code256 (with the
# code of bigsnpr, 0, 1, 2 and NA), saved next to its backing file; this
# replaces bigsnpr::snp_readBed(), which also reads the .bim and .fam files
# (here they are read by loci_from_bim() and indiv_meta_from_fam()). The
# genotypes are decoded in parallel with n_cores cores (read_bed_to_fbm(),
# in src/read_bed.cpp). Returns the path of the .rds file of the FBM.
fbm_read_bed <- function(bedfile,
                         n_indiv,
                         n_snp,
                         backingfile = bigsnpr::sub_bed(bedfile),
                         n_cores = 1) {
  backingfile <- path.expand(backingfile)
  bedfile <- path.expand(bedfile)
  # check that file exists
  if (file.exists(paste0(backingfile, ".bk"))) {
    stop("the backing file ", paste0(backingfile, ".bk"), " already exists")
  }
  if (file.exists(paste0(backingfile, ".rds"))) {
    stop("the backing file ", paste0(backingfile, ".rds"), " already exists")
  }
  # check the .bed file before creating the backing file (which can be very
  # large), so that a bad file leaves nothing behind
  check_bed_file(bedfile, n_indiv = n_indiv, n_snp = n_snp)
  big_geno <- bigstatsr::FBM.code256(
    nrow = n_indiv, ncol = n_snp,
    code = bigsnpr::CODE_012, backingfile = backingfile, init = NULL,
    create_bk = TRUE
  )
  read_bed_to_fbm(bedfile, big_geno, getCode(), n_cores)
  rds <- bigstatsr::sub_bk(big_geno$backingfile, ".rds")
  saveRDS(big_geno, rds)
  rds
}

# Check that a PLINK .bed file starts with the magic number of a SNP-major
# file, and that its size is exactly that of the magic number (3 bytes) and of
# n_snp loci of ceiling(n_indiv / 4) bytes each, as PLINK requires. A larger
# file is an error too: the .bim or .fam file does not match it (e.g. it was
# truncated).
check_bed_file <- function(bedfile, n_indiv, n_snp) {
  magic <- readBin(bedfile, what = "raw", n = 3)
  if (length(magic) < 3 || magic[1] != as.raw(0x6c) ||
        magic[2] != as.raw(0x1b)) {
    stop(
      "the file ", bedfile, " is not a PLINK .bed file (wrong magic number)"
    )
  }
  if (magic[3] != as.raw(0x01)) {
    stop(
      "the .bed file ", bedfile, " is in individual-major mode, which is ",
      "not supported; convert it to SNP-major mode with PLINK (e.g. plink ",
      "--bfile <prefix> --make-bed --out <new prefix>)"
    )
  }
  # (sizes as doubles, which are exact up to 2^53 bytes)
  expected_size <- 3 + as.numeric(n_snp) * ceiling(n_indiv / 4)
  bed_size <- file.size(bedfile)
  if (bed_size < expected_size) {
    stop(
      "the .bed file ", bedfile, " is smaller than expected for the number ",
      "of individuals in the .fam file (", n_indiv, ") and of loci in the ",
      ".bim file (", n_snp, "); it might be truncated or corrupted"
    )
  }
  if (bed_size > expected_size) {
    stop(
      "the .bed file ", bedfile, " is larger than expected for the number ",
      "of individuals in the .fam file (", n_indiv, ") and of loci in the ",
      ".bim file (", n_snp, "); the .bim or .fam file might be truncated, ",
      "or the files might not belong together"
    )
  }
  invisible(TRUE)
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
#' @param names_as_int logical, if TRUE the names of the loci are their
#'   index in the bim file (an integer), and the names in the file are not
#'   read.
#' @returns A tibble with columns: `big_index`, `name`, `chromosome`,
#'   `position`, `gentic_dist`, `allele_ref`, `allele_alt`
#' @keywords internal
#' @noRd

loci_from_bim <- function(bim,
                          valid_alleles = c("A", "T", "C", "G"),
                          missing_alleles = c("0", "."),
                          allow_duplicates = FALSE,
                          names_as_int = FALSE) {
  if (is.character(bim)) {
    if (names_as_int) {
      # the names are not read (they can take several GB for tens of
      # millions of loci)
      bim <- read_loci_file(bim, character_cols = c(5, 6), drop = 2)
      names(bim)[1:5] <- c(
        "chromosome", "genetic.dist", "physical.pos", "allele1", "allele2"
      )
    } else {
      bim <- read_loci_file(bim, character_cols = c(2, 5, 6))
      names(bim)[1:6] <- c(
        "chromosome", "marker.ID", "genetic.dist",
        "physical.pos", "allele1", "allele2"
      )
    }
  }
  if (is.data.frame(bim) && names_as_int) {
    # the names of the loci are their index (as for vcf files)
    bim$marker.ID <- seq_len(nrow(bim))
  }
  if (!is.data.frame(bim) || ncol(bim) < 6) {
    stop(
      "bim must be a data.frame with at least 6 columns or a path",
      "to a .bim file"
    )
  }
  loci <- tibble::tibble(
    name = if (names_as_int) bim$marker.ID else as.character(bim$marker.ID),
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

#' Read a whitespace delimited file with one row per locus
#'
#' This function reads files such as PLINK .bim or *packedancestry* .snp files
#' with [data.table::fread()], which is much faster than [utils::read.table()]
#' and uses much less memory on large files (read.table() creates a string for
#' every field before converting the columns, which costs several GB for tens
#' of millions of loci).
#' @param file the path to the file
#' @param character_cols the indices of the columns to read as characters
#'   (locus names and alleles, so that e.g. an allele "T" is not read as
#'   TRUE, nor a name "001" as the integer 1); the type of the other columns
#'   is guessed, as done by read.table() (e.g. chromosomes are integers if
#'   all values are integers, so that their factor levels are sorted
#'   numerically)
#' @param drop the indices of the columns not to read (`NULL` to read all)
#' @returns a data.frame with columns V1, V2, ... (named after their position
#'   in the file)
#' @keywords internal
#' @noRd
read_loci_file <- function(file, character_cols, drop = NULL) {
  data.table::fread(
    file,
    header = FALSE,
    drop = drop,
    colClasses = list(character = character_cols),
    # read large integers (e.g. positions) as doubles, not integer64
    integer64 = "double",
    data.table = FALSE,
    showProgress = FALSE
  )
}


#' Generate a tibble of indiv_meta from a fam file
#'
#' This function reads a PLINK .fam file and extracts the relevant information
#' to create an indiv_meta tibble compatible with a `gen_tibble` object.
#' @param fam the path to a fam file or a data.frame containing the contents of
#' a fam file.
#' @returns A data.frame with columns: `id`, `population`, `paternal_ID`,
#' `maternal_ID`, `sex`, `phenotype` (a factor for case/control phenotypes,
#' or numeric for quantitative ones); `population` is only present if the
#' family ID differs from the individual ID for some individuals, and the
#' other columns if they have some non-missing values
#' @keywords internal
#' @noRd

indiv_meta_from_fam <- function(fam) {
  if (is.character(fam)) {
    # read the family, individual and parental IDs as characters, so that e.g.
    # "001" is not turned into "1" (the type of sex and phenotype is guessed)
    fam <- utils::read.table(
      fam,
      stringsAsFactors = FALSE,
      colClasses = c(rep("character", 4), NA, NA)
    )
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
      indiv_meta$sex <- dplyr::recode_values(
        fam$sex,
        1 ~ "male",
        2 ~ "female",
        default = NA,
        ptype = factor(levels = c("female", "male"))
      )
    }
  }

  # as in PLINK, the phenotype is case/control (1 = control, 2 = case, and 0
  # or -9 = missing) unless it has other values, in which case it is
  # quantitative (with -9 as missing)
  if (is.numeric(fam$affection) && !all(fam$affection %in% c(0, -9))) {
    if (all(fam$affection %in% c(0, 1, 2, -9, NA))) {
      indiv_meta$phenotype <- dplyr::recode_values(
        fam$affection,
        1 ~ "control",
        2 ~ "case",
        default = NA,
        ptype = factor(levels = c("control", "case"))
      )
    } else {
      indiv_meta$phenotype <- as.numeric(fam$affection)
      indiv_meta$phenotype[indiv_meta$phenotype %in% -9] <- NA
    }
  }

  indiv_meta <- validate_indiv_meta(as.data.frame(indiv_meta))

  return(indiv_meta)
}
