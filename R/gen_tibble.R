#' Constructor for a `gen_tibble`
#'
#' A `gen_tibble` stores genotypes for individuals in a tidy format. DESCRIBE
#' here the format
#'
#' - *VCF* files: the fast `cpp` parser is used by default. Both `cpp` and
#' `vcfR` parsers attempt to establish ploidy from the first variant; if that
#' variant is found in a sex chromosome (or mtDNA), the parser will fail with
#' 'Error: a genotype has more than max_ploidy alleles...'. To successful import
#' such a *VCF*, change the order of variants so that the first chromosome is an
#' autosome using a tool such as `vcftools`. Currently, only biallelic SNPs are
#' supported. If haploid variants (e.g. sex chromosomes) are included in the
#' *VCF*, they are not transformed into homozygous calls. Instead, reference
#' alleles will be coded as 0 and alternative alleles will be coded as 1.
#'
#' - *packedancestry* files: When loading *packedancestry* files,
#' missing alleles will be converted from 'X' to NA
#'
#' @note Helper functions for accessing `gen_tibble` object attributes and
#' checking gen_tibble ploidy can be found in gt_helper_functions.R
#'
#' @param x can be:
#' - a string giving the path to a PLINK BED or PED file. The associated
#'   BIM and FAM files for the BED, or MAP for PED are expected to be in the
#'   same directory and have the same file name.
#' - a string giving the path to a RDS file storing a `bigSNP` object from
#'   the `bigsnpr` package (usually created with [bigsnpr::snp_readBed()])
#' - a string giving the path to a vcf file. Only biallelic SNPs will be
#'   considered.
#' - a string giving the path to a *packedancestry* .geno file. The associated
#'   .ind and .snp files are expected to be in the same directory and share the
#'   same file name prefix.
#' - a genotype matrix of dosages (0, 1, 2, NA) giving the dosage of the
#'   alternate allele.
#' @param indiv_meta a list, data.frame or tibble with compulsory columns 'id'
#'   and 'population', plus any additional metadata of interest. This is only
#'   used if `x` is a genotype matrix. Otherwise this information is extracted
#'   directly from the files.
#' @param loci a data.frame or tibble, with compulsory columns 'name',
#'   'chromosome', and 'position','genetic_dist', 'allele_ref' and 'allele_alt'.
#'   This is only used if `x` is a genotype matrix. Otherwise this information
#'   is extracted directly from the files.
#' @param chunk_size the number of loci or individuals (depending on the format)
#'   processed at a time (currently used if `x` is a vcf with parser "vcfR")
#' @param ... if `x` is the name of a vcf file, additional arguments passed to
#'   [vcfR::read.vcfR()]. Otherwise, unused.
#' @param parser the name of the parser used for *VCF*, either "cpp" to use a
#'   fast C++ parser (the default), or "vcfR" to use the R package `vcfR`. The
#'   latter is slower but more robust; if "cpp" gives an error, try using "vcfR"
#'   in case your
#'   *VCF* has an unusual structure.
#' @param n_cores the number of cores to use for parallel processing
#' @param valid_alleles a vector of valid allele values; it defaults to 'A','T',
#'   'C' and 'G'.
#' @param missing_alleles a vector of values in the BIM file/loci dataframe that
#'   indicate a missing value for the allele value (e.g. when we have a
#'   monomorphic locus with only one allele). It defaults to '0' and '.' (the
#'   same as PLINK 1.9).
#' @param backingfile the path, including the file name without extension, for
#'   backing files used to store the data (they will be given a .bk and .RDS
#'   automatically). This is not needed if `x` is already an .RDS file. If `x`
#'   is a .BED or a *VCF* file and `backingfile` is left NULL, the backing file
#'   will be saved in the same directory as the bed or vcf file, using the same
#'   file name but with a different file type (.bk rather than .bed or .vcf). If
#'   `x` is a genotype matrix and `backingfile` is NULL, then a temporary file
#' will be created (but note that R will delete it at the end of the session!)
#' @param allow_duplicates logical. If TRUE, the tibble will allow duplicated
#'   loci (those with genomic coordinate (chromosome + position) or locus name
#'   appearing more than once). If FALSE, an error will be thrown if duplicated
#'   loci are found. These validations run before backing
#'   files are saved. Default is FALSE.
#' @param quiet provide information on the files used to store the data
#' @returns an object of the class `gen_tbl`.
#' @rdname gen_tibble
#' @export
#' @examples
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' # Create a gen_tibble from a .bed file
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Create a gen_tibble from a .vcf file
#' vcf_path <-
#'   system.file("extdata", "anolis",
#'     "punctatus_t70_s10_n46_filtered.recode.vcf.gz",
#'     package = "tidypopgen"
#'   )
#' gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile("anolis_"))
#'
#' # Create a gen_tibble from a matrix of genotypes:
#' test_indiv_meta <- data.frame(
#'   id = c("a", "b", "c"),
#'   population = c("pop1", "pop1", "pop2")
#' )
#' test_genotypes <- rbind(
#'   c(1, 1, 0, 1, 1, 0),
#'   c(2, 1, 0, 0, 0, 0),
#'   c(2, 2, 0, 0, 1, 1)
#' )
#' test_loci <- data.frame(
#'   name = paste0("rs", 1:6),
#'   chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
#'   position = as.integer(c(3, 5, 65, 343, 23, 456)),
#'   genetic_dist = as.double(rep(0, 6)),
#'   allele_ref = c("A", "T", "C", "G", "C", "T"),
#'   allele_alt = c("T", "C", NA, "C", "G", "A")
#' )
#'
#' gen_tibble(
#'   x = test_genotypes,
#'   loci = test_loci,
#'   indiv_meta = test_indiv_meta,
#'   valid_alleles = c("A", "T", "C", "G"),
#'   quiet = TRUE
#' )
#'
gen_tibble <-
  function(x,
           ...,
           valid_alleles = c("A", "T", "C", "G"),
           missing_alleles = c("0", "."),
           backingfile = NULL,
           allow_duplicates = FALSE,
           quiet = FALSE) {
    UseMethod("gen_tibble", x)
  }

###############################################################################
# character method for files
###############################################################################
#' @export
#' @rdname gen_tibble
gen_tibble.character <-
  function(x,
           ...,
           parser = c("cpp", "vcfR"),
           n_cores = 1,
           chunk_size = NULL,
           valid_alleles = c("A", "T", "C", "G"),
           missing_alleles = c("0", "."),
           backingfile = NULL,
           allow_duplicates = FALSE,
           quiet = FALSE) {
    # parser for vcf
    parser <- match.arg(parser)
    # check that valid alleles does not contain zero
    if ("0" %in% valid_alleles) {
      stop(paste(
        "'0' cannot be a valid allele",
        "(it is the default missing allele value!)"
      ))
    }

    if (is.null(backingfile)) {
      backingfile <- change_duplicated_file_name(x)
    } else if (!is.null(backingfile)) {
      backingfile <- change_duplicated_file_name(backingfile)
    }


    # test that the string is not empty (or it will break the code later on)
    if (nchar(x) == 0) {
      stop("x should not be an empty string")
    }

    if ((tolower(file_ext(x)) == "bed") || (tolower(file_ext(x)) == "rds")) {
      rlang::check_dots_empty()
      x_gt <- gen_tibble_bed_rds(
        x = x,
        ...,
        valid_alleles = valid_alleles,
        missing_alleles = missing_alleles,
        backingfile = backingfile,
        allow_duplicates = allow_duplicates,
        quiet = quiet
      )
    } else if (
      (tolower(file_ext(x)) == "vcf") || (tolower(file_ext(x)) == "gz")
    ) {
      # note that gen_tibble_vcf generates the files for a bigSNP object, which
      # is then passed back to gen_tibble_bed_rds to create the gen_tibble
      # so, the object returned by gen_tibble_vcf is already the
      # final gen_tibble
      x_gt <- gen_tibble_vcf(
        x = x,
        ...,
        parser = parser,
        chunk_size = chunk_size,
        n_cores = n_cores,
        valid_alleles = valid_alleles,
        missing_alleles = missing_alleles,
        backingfile = backingfile,
        allow_duplicates = allow_duplicates,
        quiet = quiet
      )
    } else if (tolower(file_ext(x)) == "ped") {
      x_gt <- gen_tibble_ped(
        x = x,
        ...,
        valid_alleles = valid_alleles,
        missing_alleles = missing_alleles,
        backingfile = backingfile,
        allow_duplicates = allow_duplicates,
        quiet = quiet
      )
    } else if (tolower(file_ext(x)) == "geno") {
      x_gt <- gen_tibble_packedancestry(
        x = x,
        ...,
        valid_alleles = valid_alleles,
        missing_alleles = missing_alleles,
        backingfile = backingfile,
        allow_duplicates = allow_duplicates,
        quiet = quiet,
        chunk_size = chunk_size
      )
    } else {
      stop(paste(
        "x should be a valid file path pointing to a either a PLINK .bed",
        "or .ped file, a bigSNP .rds file or a VCF .vcf or",
        ".vcf.gz file"
      ))
    }

    # check alleles
    loci <- show_loci(x_gt)
    doubles <- which(is.na(loci$allele_ref) & is.na(loci$allele_alt))

    if (length(doubles) > 0) {
      check_missing <- x_gt %>%
        select_loci(all_of(doubles)) %>%
        loci_missingness()
      if (any(check_missing < 1, na.rm = TRUE)) {
        # clean up files if we are stopping
        files <- gt_get_file_names(x_gt)
        if (file.exists(files[1])) file.remove(files[1])
        if (file.exists(files[2])) file.remove(files[2])
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

    gt_save_light(x_gt, quiet = quiet) # nolint
    return(x_gt)
  }


###############################################################################
# matrix method to provide data directly from R
###############################################################################
#' @param ploidy the ploidy of the samples (either a single value, or
#' a vector of values for mixed ploidy). Only used if creating
#' a gen_tibble from a matrix of data; otherwise, ploidy is determined
#' automatically from the data as they are read.
#' @export
#' @rdname gen_tibble
gen_tibble.matrix <- function(
    x,
    indiv_meta,
    loci,
    ...,
    ploidy = 2,
    valid_alleles = c("A", "T", "C", "G"),
    missing_alleles = c("0", "."),
    backingfile = NULL,
    allow_duplicates = FALSE,
    quiet = FALSE) {
  rlang::check_dots_empty()

  # check that valid alleles does not contain zero
  if ("0" %in% valid_alleles) {
    stop(paste(
      "'0' cannot be a valid allele",
      "(it is the default missing allele value!)"
    ))
  }

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
    max_ploidy <- max(ploidy)
  } else {
    if ((ploidy != 0) && (ploidy != -2)) {
      fbm_ploidy <- rep(ploidy, nrow(x))
    } else {
      stop(
        "'ploidy' 0 (mixed ploidy) or -2 (haplodiploids) ",
        "require a vector of individual ploidies"
      )
    }
    max_ploidy <- ploidy
  }

  loci <- validate_loci(loci,
    check_alphabet = TRUE,
    check_duplicates = TRUE,
    allow_duplicates = allow_duplicates,
    harmonise_loci = TRUE,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles
  )
  indiv_meta <- validate_indiv_meta(indiv_meta)

  # validate x (the genotypes)
  # check that x (the genotypes) is numeric matrix
  if (inherits(x, "data.frame")) {
    x <- as.matrix(x)
  }
  if (!(is.matrix(x) && is.numeric(x))) {
    stop("'x' is not a numeric matrix of integers")
  }

  # check dimensions
  if (ncol(x) != nrow(loci)) {
    stop(paste(
      "there is a mismatch between the number of loci in the",
      "genotype table x and in the loci table"
    ))
  }
  if (nrow(x) != nrow(indiv_meta)) {
    stop(paste(
      "there is a mismatch between the number of individuals in the",
      "genotype table x and in the indiv_meta table"
    ))
  }

  if (!is.null(backingfile)) {
    backingfile <- change_duplicated_file_name(backingfile)
  }

  fbm_obj <- gt_write_fbm_from_dfs(
    genotypes = x,
    backingfile = backingfile,
    max_ploidy = max_ploidy
  )

  fbm_path <- bigstatsr::sub_bk(fbm_obj$backingfile, ".rds")

  indiv_id <- indiv_meta$id

  indiv_meta <- as.list(indiv_meta)

  indiv_meta$genotypes <- new_vctrs_bigsnp(
    fbm_obj,
    fbm_file = fbm_path,
    loci = loci,
    indiv_id = indiv_id,
    ploidy = max_ploidy,
    fbm_ploidy = fbm_ploidy
  )

  new_gen_tbl <- tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )

  # check alleles
  loci <- show_loci(new_gen_tbl)
  doubles <- which(is.na(loci$allele_ref) & is.na(loci$allele_alt))

  if (length(doubles) > 0) {
    check_missing <- new_gen_tbl %>%
      select_loci(all_of(doubles)) %>%
      loci_missingness()
    if (any(check_missing < 1, na.rm = TRUE)) {
      files <- gt_get_file_names(new_gen_tbl)
      if (file.exists(files[1])) file.remove(files[1])
      if (file.exists(files[2])) file.remove(files[2])
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

  gt_save(new_gen_tbl, quiet = quiet) # nolint
  return(new_gen_tbl)
}


################################################################################
## vctrs_bigSNP class to store the genotype info in a gen_tibble
################################################################################

#' create a vctrs_bigSNP
#' @param fbm_obj the FBM object (bigstatsr::FBM.code256)
#' @param fbm_file the rds file associated with the FBM object
#' @param loci a tibble of loci (needs to be validated first with
#'   `validate_loci`)
#' @param indiv_id a vector of individual ids (from indiv_meta)
#' @param ploidy the ploidy of the samples, a single value.
#' @param fbm_ploidy a vector of ploidies for each individual in the fbm object
#'   (note that this is for the full FBM, not just the tibble)
#' @returns a vctrs_bigSNP object
#' @keywords internal
#' @noRd
new_vctrs_bigsnp <- function(fbm_obj, fbm_file, loci, indiv_id, ploidy = 2,
                             fbm_ploidy = NULL) {
  # check that indiv_id is the same length as the nrow of fmb_obj
  if (length(indiv_id) != nrow(fbm_obj)) {
    stop(paste(
      "'indiv_id' should be the same length as the number of rows",
      "in the fbm_obj"
    ))
  }
  # check that nrow(loci) is ncol(fbm_obj)
  if (nrow(loci) != ncol(fbm_obj)) {
    stop(paste(
      "'loci' should have the same number of rows as the number of",
      "columns in the fbm_obj"
    ))
  }

  if (ploidy != -2) {
    if (length(unique(fbm_ploidy)) > 1) {
      max_ploidy <- 0
    } else {
      max_ploidy <- max(fbm_ploidy)
    }
  } else {
    max_ploidy <- -2
  }

  # add the big_index column
  loci <- loci %>% dplyr::mutate(big_index = dplyr::row_number(), .before = 1)

  new_vectr <- vctrs::new_vctr(
    seq_along(indiv_id),
    fbm = fbm_obj,
    fbm_file = fbm_file,
    # TODO make sure this does not take too long
    fbm_md5sum = tools::md5sum(fbm_file),
    loci = loci,
    names = indiv_id,
    ploidy = max_ploidy,
    class = "vctrs_bigSNP"
  )
  attr(new_vectr, "fbm_ploidy") <- fbm_ploidy
  new_vectr
}

#' @export
summary.vctrs_bigSNP <- function(object, ...) {
  summary(rep("bigSNP-genotypes", length(object)))
}


################################################################################

#' Test if a tibble is really `gen_tibble`
#'
#' Some `dplyr` operations strip the subclass from the tibble. This function
#' is used to check if the tibble is, in reality, still of class `gen_tbl`
#' @param .x the tibble
#' @returns the gen_tibble, invisibly
#' @keywords internal
#' @noRd
stopifnot_gen_tibble <- function(.x) {
  if (!"genotypes" %in% names(.x)) {
    stop("not a gen_tibble, 'genotypes' column is missing")
  }
  if (!inherits(.x$genotypes, "vctrs_bigSNP")) {
    stop("not a gen_tibble, the genotypes column is not of class vctrs_bigSNP")
  }
  return(invisible(.x))
}

# print method
#' @export
tbl_sum.gen_tbl <- function(x, ...) {
  c(
    "A gen_tibble" = paste(count_loci(x), " loci"),
    NextMethod()
  )
}


# function to check the allele alphabet
check_allele_alphabet <- function(
    x,
    valid_alleles = c("A", "T", "C", "G"),
    missing_alleles = c("0", ".")) {
  if (
    any(
      !x$allele_ref %in% c(valid_alleles, missing_alleles, NA),
      !x$allele_alt %in% c(valid_alleles, missing_alleles, NA)
    )
  ) {
    stop(
      "valid alleles are ",
      paste(c(valid_alleles, missing_alleles), collapse = " "),
      " but ",
      paste(
        unique(c(
          x$allele_ref,
          x$allele_alt
        )),
        collapse = " "
      ),
      " were found."
    )
  }
}

# set all missing values to NA
# loci_info is usually from show_loci()
harmonise_missing_values <- function(loci_info, missing_alleles = c("0", ".")) {
  # 0 is always considered as a missing value
  if (!"0" %in% missing_alleles) {
    missing_alleles <- c(missing_alleles, "0")
  }
  loci_info$allele_ref[loci_info$allele_ref %in% missing_alleles] <- NA
  loci_info$allele_alt[loci_info$allele_alt %in% missing_alleles] <- NA
  return(loci_info)
}


# check for existing .bk files
change_duplicated_file_name <- function(file) {
  file <- tools::file_path_sans_ext(file)

  bk <- paste0(file, ".bk")
  rds <- paste0(file, ".rds")

  if (file.exists(bk) || file.exists(rds)) {
    version <- 2
    # extract the base name and version number
    base_name_pattern <- "^(.*)_v(\\d+)$"
    # check for any matches
    matches <- regmatches(
      basename(file),
      regexec(base_name_pattern, basename(file))
    )

    if (length(matches[[1]]) > 0) {
      # Extract base name without version and current version number
      base_name <- matches[[1]][2] # Part before "_v<number>"
      # extract current version number
      current_version <- as.numeric(matches[[1]][3]) # nolint
    } else {
      base_name <- basename(file) # Use the full name if there's no "_v" suffix
      current_version <- 1
    }

    version_pattern <- paste0(base_name, "_v(\\d+)\\.bk$")
    # read files to check for existing versions
    existing_files <- list.files(
      dirname(bk),
      pattern = paste0(
        "^",
        base_name,
        "_v\\d+\\.bk$"
      )
    )

    if (length(existing_files) > 0) {
      versions <- sub(version_pattern, "\\1", existing_files)
      versions <- as.numeric(versions)
      if (!any(is.na(versions))) {
        version <- max(versions) + 1 # add 1 to the version number
      }
    }
    # create new file path
    new_file <- file.path(dirname(file), paste0(base_name, "_v", version))

    return(new_file)
  }
  return(file) # nolint
}

# adding a chr_int column
cast_chromosome_to_int <- function(chromosome) {
  # if chromosome is an factor, then cast it to character
  if (is.factor(chromosome)) {
    chromosome <- as.character(chromosome)
  }
  # if chromosome is a character, then cast it to integer
  if (is.character(chromosome)) {
    # attempt to strip chr from the chromosome
    chromosome <- gsub(
      "^(chromosome_|chr_|chromosome|chr)",
      "",
      chromosome,
      ignore.case = TRUE
    )
    chromosome <- tryCatch(as.integer(chromosome), warning = function(e) {
      as.integer(as.factor(chromosome))
    })
  }
  if (is.numeric(chromosome)) {
    chromosome <- as.integer(chromosome)
  }
  if (is.integer(chromosome)) {
    return(chromosome) # nolint
  } else {
    stop("Chromosome column should be integer, character, or factor")
  }
}
