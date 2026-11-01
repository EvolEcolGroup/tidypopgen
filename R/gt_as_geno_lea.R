#' Convert a `gentibble` to a .geno file for sNMF from the LEA package
#'
#' This function writes a .geno file from a [`gen_tibble`]. Unless a file path
#' is given, a file with suffix .geno is written in the same location as the
#' .rds and .bk files that underpin the [`gen_tibble`].
#'
#' The genotypes are read and written directly from the file backed matrix
#' (one locus at a time), so this function also works on datasets that do not
#' fit in memory.
#'
#' @note the `.geno` format only supports haploid or diploid data, so this
#'   function only works on diploid (or pseudohaploid) `gen_tibble` objects.
#'   Pseudohaploid data is simply treated as regular diploid genotype counts
#'   (dosage 0/2), the same convention used by other tools (e.g. PLINK).
#' @param x a [`gen_tibble`]
#' @param file the .geno filename with a path, or NULL (the default) to use the
#'   location of the backing files.
#' @returns the path of the .geno file
#' @export
#' @seealso [LEA::geno()]
#' @examplesIf rlang::is_installed("LEA")
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Write a geno file
#' gt_as_geno_lea(example_gt, file = paste0(tempfile(), "_example.geno"))
# Write a geno file (as used by LEA)
gt_as_geno_lea <- function(x, file = NULL) {
  stopifnot_dip_pseudo(x)
  if (is.null(file)) {
    file <- sub_rds(gt_get_file_names(x)[1], ".geno")
  }
  if (file_ext(file) != "geno") {
    file <- paste0(file, ".geno")
  }
  # create a lookup table with the string for each raw value of the FBM, from
  # its decoded dosages (missing values are coded as 9)
  genotypes_lut <- .gt_get_fbm(x)$code256
  genotypes_lut[is.na(genotypes_lut)] <- 9
  genotypes_lut <- as.character(genotypes_lut)
  # create an empty file (overwriting any existing file), and append the
  # genotypes to it in C++, one line per locus (with no separator between
  # individuals)
  file.create(file)
  gt_write_text_by_locus(
    BM = .gt_get_fbm(x),
    rowInd = .gt_fbm_rows(x),
    colInd = .gt_fbm_cols(x),
    lut = genotypes_lut,
    prefix = character(0),
    sep = "",
    file = path.expand(file)
  )
  # return the file name
  return(file)
}
