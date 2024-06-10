#' Save a gen_tibble
#'
#' Save the tibble (and update the backing files). A full explanation of files (names,extensions, etc.)
#'
#' @param x a [`gen_tibble`]
#' @param file_name the file name, including the full path. If it does not end
#' with '.gt', the extension will be added.
#' @param quiet boolean to suppress information about hte files
#' @returns the file name and path of the .gt file, together with the .rds and .bk files
#' @export

gt_save <-function(x, file_name = NULL, quiet = FALSE) {
  if (!inherits(x, "gen_tbl")){
    stop("x should be a gen_tibble")
  }

  if (is.null(file_name)){
    file_name <-  bigstatsr::sub_bk(gt_get_file_names(x)[2],".gt")
  }
  if (file_ext(file_name)!="gt"){
    file_name <- paste0(file_name,".gt")
  }
  # we update the bigsnp object
  bigsnpr::snp_save(attr(x$genotypes,"bigsnp"))
  # and now save our gen_tibble
  saveRDS(x,file_name)
  if (!quiet){
    message("\ngen_tibble saved to ", file_name)
    message("using bigSNP file: ", gt_get_file_names(x)[1])
    message("with backing file: ", gt_get_file_names(x)[2])
    message("make sure that you do NOT delete those files!")
    message("to reload the gen_tibble in another session, use gt_load('",file_name,"')")
  }
  return(c(file_name,gt_get_file_names(x)))
}


sub_gt <- function (path, replacement = "", stop_if_not_ext = TRUE)
{
  pattern <- "\\.gt$"
  if (!grepl(pattern, path))
    stop("Path '%s' must have 'gt' extension.", path)
  if (stop_if_not_ext && (nchar(replacement) > 0) && (substr(replacement,
                                                             1, 1) != "."))
    stop("Replacement must be an extension starting with '.' if provided.")
  sub(pattern, replacement, path)
}

#' Get the names of files storing the genotypes of a `gen_tibble`
#'
#' A function to return the names of the files used to store data in a `gen_tibble`.
#' Specifically, this returns the .rds file storing the big
#' @param x a [`gen_tibble`]
#' @returns a character vector with the names and paths of the two files
#' @export
gt_get_file_names <- function(x){
  return(c(attr(x$genotypes,"bigsnp")$genotypes$rds,
    attr(x$genotypes,"bigsnp")$genotypes$backingfile))
}
