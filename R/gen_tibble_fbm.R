#' Test that a loci table is valid
#' 
#' This function checks that a loci tibble has the required columns and that they
#' are of the correct type.
#' @param loci A tibble of loci
#' @returns the validated loci table
#' @keywords internal
#' @noRd

validate_loci <- function(loci) {
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
  # if chromosome is numeric, turn it into character
  if (is.numeric(loci$chromosome)) {
    loci$chromosome <- as.character(loci$chromosome)
  }
  # check that chromosome is a factor or character
  if (!is.factor(loci$chromosome) && !is.character(loci$chromosome)) {
    stop("loci$chromosome must be a factor or character")
  }
  if (!is.integer(loci$position) && !is.numeric(loci$position)) {
    stop("loci$position must be an integer")
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
  return(loci)
}


#' Validate indiv_meta
#' 
#' This function checks that an indiv_meta tibble has the required columns and that they
#' are of the correct type (just the id column).
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
    ploidy = 2) {
  if (is.null(backingfile)) {
    backingfile <- tempfile()
  }
  # set up code (accounting for ploidy)
  code256 <- rep(NA_real_, 256)
  if (length(ploidy) > 1) {
    # check that there are no missing values in ploidy vector
    if (any(is.na(ploidy))) {
      stop("'ploidy' can not contain NAs")
    }
    max_ploidy <- max(ploidy)
  } else {
    max_ploidy <- ploidy
  }
  
  if (is.na(max_ploidy)) {
    stop("'ploidy' can not contain NAs")
  }
  
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
