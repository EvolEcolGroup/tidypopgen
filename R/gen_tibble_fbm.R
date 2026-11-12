#' Test that a loci table is valid
#'
#' This function checks that a loci tibble has the required columns and that
#' they are of the correct type.
#' @param loci A tibble of loci
#' @param check_alphabet whether to check that alleles are in the valid_alleles
#'   list. Default is FALSE.
#' @param harmonise_loci whether to harmonise missing alleles in the loci table
#'   using the missing_alleles list. Default is FALSE.
#' @param check_duplicates whether to check for duplicated loci (same chromosome
#'   and position) or duplicated locus names. Default is FALSE.
#' @param allow_duplicates whether to allow duplicated loci (same chromosome and
#'   position) or duplicated locus names. Default is FALSE.
#' @param valid_alleles a character vector of valid alleles. Default is c("A",
#'   "T", "C", "G").
#' @param missing_alleles a character vector of alleles to be treated as
#'   missing. Default is c("0", ".").
#' @returns the validated loci table
#' @keywords internal
#' @noRd

validate_loci <- function(loci,
                          check_alphabet = FALSE,
                          harmonise_loci = FALSE,
                          check_duplicates = FALSE,
                          allow_duplicates = FALSE,
                          valid_alleles = c("A", "T", "C", "G"),
                          missing_alleles = c("0", "."),
                          chromosomes_as_int = FALSE) {
  required_cols <- c(
    "name", "chromosome", "position",
    "genetic_dist", "allele_ref", "allele_alt"
  )
  if (!is.data.frame(loci)) {
    stop("loci must be a data.frame or a tibble")
  }
  if (!all(required_cols %in% colnames(loci))) {
    stop(paste0(
      "loci must have the following columns: ",
      paste(required_cols, collapse = ", ")
    ))
  }
  # loci names should be a character or an integer
  if (!is.character(loci$name) && !is.integer(loci$name)) {
    stop("loci$name must be a character or an integer")
  }
  # chromosomes are a factor, unless we keep their integer codes
  # (chromosomes_as_int: with millions of chromosomes, a factor would create
  # a string for each of them)
  if (chromosomes_as_int) {
    if (!is.numeric(loci$chromosome)) {
      stop("loci$chromosome must be integer codes")
    }
    if (!is.integer(loci$chromosome)) {
      loci$chromosome <- as.integer(loci$chromosome)
    }
  } else if (!is.factor(loci$chromosome)) {
    loci$chromosome <- as.factor(loci$chromosome)
  }
  if (!is.integer(loci$position) && !is.numeric(loci$position)) {
    stop("loci$position must be integer-like (integer or numeric)")
  }
  if (!is.numeric(loci$genetic_dist)) {
    stop("loci$genetic_dist must be numeric")
  }
  if (!is.character(loci$allele_ref)) {
    stop("loci$allele_ref must be a character")
  }
  if (!is.character(loci$allele_alt)) {
    stop("loci$allele_alt must be a character")
  }
  if (check_alphabet == TRUE) {
    check_allele_alphabet(loci,
      valid_alleles = valid_alleles,
      missing_alleles = missing_alleles
    )
  }
  if (harmonise_loci == TRUE) {
    loci <- harmonise_missing_values(
      loci,
      missing_alleles
    )
  }
  if (check_duplicates == TRUE) {
    # check for duplicated positions within chromosomes, i.e. duplicated
    # pairs of chromosome and position. When positions are sorted within
    # each chromosome (as usual), duplicates are equal consecutive positions,
    # which a single pass in C++ finds without hashing every pair (each
    # chromosome is a group of its own); otherwise (or with missing
    # positions) we hash the pairs with vctrs
    chromosome_groups <- chromosome_group_codes(loci$chromosome)
    position_order <- unsorted_by_chromosome_cpp(
      loci$position,
      chromosome_groups$codes,
      seq_len(chromosome_groups$n)
    )
    if (isFALSE(position_order[["unsorted"]])) {
      has_dup_pos <- position_order[["ties"]]
    } else {
      has_dup_pos <- vctrs::vec_duplicate_any(
        vctrs::data_frame(
          chromosome = loci$chromosome,
          position = loci$position
        )
      )
    }
    has_dup_names <- any_duplicated_names(loci$name) > 0

    if (!allow_duplicates) {
      if (has_dup_pos) {
        stop(paste0(
          "Your data contain duplicated loci. ",
          "Remove them or set allow_duplicates = TRUE."
        ))
      }
      if (has_dup_names) {
        stop(paste0(
          "Your data contain duplicated locus names. ",
          "Remove them or set allow_duplicates = TRUE."
        ))
      }
    } else {
      if (has_dup_pos) {
        warning(paste0(
          "You have allowed duplicated loci in your data. ",
          "Your data contain duplicated loci. ",
          "Use find_duplicated_loci(my_tibble) to select and remove them."
        ))
      }
      if (has_dup_names) {
        warning(paste0(
          "You have allowed duplicated loci in your data. ",
          "Your data contain duplicated locus names. ",
          "Use anyDuplicated(loci_names(my_tibble)) to select and remove them."
        ))
      }
    }
  }
  # Update chr_int column in case chromosome factor levels have changed
  loci <- tibble::as_tibble(loci)
  return(loci)
}


#' Validate indiv_meta
#'
#' This function checks that an indiv_meta tibble has the required columns and
#' that they are of the correct type (just the id column).
#' @param indiv_meta A tibble of individual metadata
#' @returns the validated indiv_meta tibble
#' @keywords internal
#' @noRd

validate_indiv_meta <- function(indiv_meta) {
  required_cols <- c("id")
  if (!is.data.frame(indiv_meta)) {
    stop("indiv_meta must be a data.frame or a tibble")
  }
  if (!all(required_cols %in% colnames(indiv_meta))) {
    stop(paste0(
      "indiv_meta must have the following columns: ",
      paste(required_cols, collapse = ", ")
    ))
  }
  # check that id is a character or an integer
  if (!is.character(indiv_meta$id) && !is.integer(indiv_meta$id)) {
    stop("indiv_meta$id must be a character or an integer")
  }
  # check that all ids are unique
  if (any(duplicated(indiv_meta$id))) {
    stop("indiv_meta$id must be unique")
  }
  return(indiv_meta)
}

#' Create an FBM object from data.frames
#'
#' This function expects the indiv_meta and loci to have the correct columns
#' @param genotypes a genotype matrix
#' @param backingfile the path, including the file name without extension, for
#'  backing files used to store the data (they will be given a .bk and .RDS
#'  automatically). If NULL, a temporary file will be created (but note that R
#'  will delete it at the end of the session!)
#' @param ploidy the ploidy of the samples (either a single value, or
#'  a vector of values for mixed ploidy).
#' @returns a FBM object
#' @keywords internal
#' @noRd
gt_write_fbm_from_dfs <- function(
  genotypes,
  backingfile = NULL,
  max_ploidy = 2
) {
  if (is.null(backingfile)) {
    backingfile <- tempfile()
  }

  # set up code (accounting for ploidy)
  code256 <- rep(NA_real_, 256)
  code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)

  # ensure max_ploidy is appropriate for the data
  dosage_check <- check_dosage_matrix_cpp(genotypes,
    ploidy = as.numeric(max_ploidy),
    check_integer = FALSE
  )
  if (dosage_check == 2) {
    stop(
      "max ploidy is set to ",
      max_ploidy,
      " but genotypes contains indviduals with greater ploidy"
    )
  }

  # equivalent to bigGeno in bigstatsr
  big_geno <- bigstatsr::FBM.code256(
    nrow = nrow(genotypes),
    ncol = ncol(genotypes),
    code = code256,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )
  # (sidecars next to the new backing file are left from older files with the
  # same name, and can not belong to the new data)
  remove_stale_sidecars(big_geno$backingfile)
  # write the dosages as bytes (missing genotypes as max_ploidy + 1) directly
  # into the FBM, without copying the matrix
  write_dosage_matrix_cpp(genotypes, big_geno, na_byte = max_ploidy + 1)

  # save it and return the path of the saved object
  rds <- bigstatsr::sub_bk(big_geno$backingfile, ".rds")
  saveRDS(big_geno, rds)
  return(big_geno)
}

# anyDuplicated() for the names of loci (a character or an integer vector);
# for character vectors, the CHARSXPs are hashed in C++, which is about 3
# times faster than anyDuplicated() and takes half of the memory (1.3 s and
# 256 MB for 40M names), unless the strings are not all in the native
# encoding (see any_duplicated_strings_cpp())
any_duplicated_names <- function(x) {
  if (is.character(x)) {
    res <- any_duplicated_strings_cpp(x)
    if (res >= 0) {
      return(res)
    }
  }
  anyDuplicated(x)
}
