#' Write a geno file for sNMF
#'
#' This function writes a geno file fom a [`gen_tibble`]. Unless a file path is given,
#' a file with suffix .geno is written in the same location as the .rds and .bk
#' file that underpin the [`gen_tibble`].
#'
#' @param x a [`gen_tibble`]
#' @param file the .geno filename with a path, or NULL (the default) to use the
#' location of the backing files.
#' @returns the path of the .geno file
#' @export
#'

# Write a geno file (as used by LEA)
gt_write_lea_geno <- function(x, file = NULL){
  if (is.null(file)){
    file <- sub_rds(gt_get_file_names(x)[1],".geno")
  }
  if (file_ext(file)!="geno"){
    file <- paste0(file,".geno")
  }
  # TODO this should be done in chunks for larger files
  geno_matrix <- t(x %>% show_genotypes())
  geno_matrix[is.na(geno_matrix)]<-9
  utils::write.table(geno_matrix ,sep="",
              file = file, row.names = FALSE,
              col.names = FALSE)
  # return the file name
  return(file)
}

