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
                          missing_alleles = c("0", ".")) {
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
  if (!is.character(loci$name)) {
    stop("loci$name must be a character")
  }
  if (!is.factor(loci$chromosome)) {
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
    # check for duplicates
    duplicated_pos <- loci %>%
      group_by(.data$chromosome) %>% # nolint
      group_map(~ .x[duplicated(.x$position) | duplicated(.x$position, fromLast = TRUE), ]$name) %>% # nolint
      unlist(use.names = FALSE)

    has_dup_pos <- length(duplicated_pos) > 0
    has_dup_names <- anyDuplicated(loci$name) > 0

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
    max_ploidy = 2) {
  if (is.null(backingfile)) {
    backingfile <- tempfile()
  }

  # set up code (accounting for ploidy)
  code256 <- rep(NA_real_, 256)
  code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)

  # ensure max_ploidy is appropriate for the data
  if (any(genotypes > max_ploidy, na.rm = TRUE)) {
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
  genotypes[is.na(genotypes)] <- max_ploidy + 1
  big_geno[] <- as.raw(genotypes)

  # save it and return the path of the saved object
  rds <- bigstatsr::sub_bk(big_geno$backingfile, ".rds")
  saveRDS(big_geno, rds)
  return(big_geno)
}
