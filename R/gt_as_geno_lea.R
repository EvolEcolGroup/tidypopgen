#' Convert a `gentibble` to a .geno file for sNMF from the LEA package
#'
#' This function writes a .geno file from a [`gen_tibble`]. Unless a file path
#' is given, a file with suffix .geno is written in the same location as the
#' .rds and .bk files that underpin the [`gen_tibble`].
#'
#' NOTE that we currently read all the data into memory to write the file, so
#' this function is not suitable for very large datasets.
#'
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
  if (is.null(file)) {
    file <- sub_rds(gt_get_file_names(x)[1], ".geno")
  }
  if (file_ext(file) != "geno") {
    file <- paste0(file, ".geno")
  }
  # TODO this should be done in chunks for larger files
  geno_matrix <- t(x %>% show_genotypes())
  geno_matrix[is.na(geno_matrix)] <- 9
  utils::write.table(
    geno_matrix,
    sep = "",
    file = file,
    row.names = FALSE,
    col.names = FALSE
  )
  # return the file name
  return(file)
}
