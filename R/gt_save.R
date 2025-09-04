#' Save a gen_tibble
#'
#' Save the tibble (and update the backing files). The `gen_tibble` object is
#' saved to a file with extension *.gt*, together with update its *.rds* and
#' *.bk* files. Note that multiple *.gt* files can be linked to the same *.rds*
#' and *.bk* files; generally, this occurs when we create multiple subsets of
#' the data. The *.gt* file then stores the information on what subset of the
#' full dataset we are interested in, whilst the *.rds* and *.bk* file store the
#' full dataset. To reload a `gen_tibble`, you can pass the name of the *.gt*
#' file with [gt_load()].
#'
#' @param x a [`gen_tibble`]
#' @param file_name the file name, including the full path. If it does not end
#'   with *.gt*, the extension will be added.
#' @param quiet boolean to suppress information about the files
#' @returns the file name and path of the *.gt* file, together with the *.rds*
#'   and *.bk* files
#' @seealso [gt_load()]
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # remove some individuals
#' example_gt <- example_gt %>% filter(id != "a")
#'
#' # save filtered gen_tibble object
#' gt_save(example_gt, file_name = paste0(tempfile(), "_example_filtered"))
#'
gt_save <- function(x, file_name = NULL, quiet = FALSE) {
  if (!inherits(x, "gen_tbl")) {
    stop("x should be a gen_tibble")
  }
  # we update the bigsnp object
  bigsnpr::snp_save(attr(x$genotypes, "bigsnp"))

  if (is.null(file_name)) {
    file_name <- bigstatsr::sub_bk(gt_get_file_names(x)[2], ".gt")
  }
  if (file_ext(file_name) != "gt") {
    file_name <- paste0(file_name, ".gt")
  }
  # and now save our gen_tibble
  saveRDS(x, file_name)
  if (!quiet) {
    message("\ngen_tibble saved to ", file_name)
    message("using bigSNP file: ", gt_get_file_names(x)[1])
    message("with backing file: ", gt_get_file_names(x)[2])
    message("make sure that you do NOT delete those files!")
    message("to reload the gen_tibble in another session, use:")
    message("gt_load('", file_name, "')")
  }
  return(c(file_name, gt_get_file_names(x)))
}

#' a light version of gt_save that does not resave the bigSNP RDS or backing
#' file, to be used internally when creating a gen_tibble (if we have just
#' created it, there is not need to resave it)
#' @param x a [`gen_tibble`]
#' @param file_name the file name, including the full path. If it does not end
#'   with *.gt*, the extension will be added.
#' @param quiet boolean to suppress information about the files
#' @returns the file name and path of the *.gt* file, together with the *.rds*
#'   and *.bk* files
#' @keywords internal
#' @noRd
gt_save_light <- function(x, file_name = NULL, quiet = FALSE) {
  if (is.null(file_name)) {
    file_name <- bigstatsr::sub_bk(gt_get_file_names(x)[2], ".gt")
  }
  if (file_ext(file_name) != "gt") {
    file_name <- paste0(file_name, ".gt")
  }
  # and now save our gen_tibble
  saveRDS(x, file_name)
  if (!quiet) {
    message("\ngen_tibble saved to ", file_name)
    message("using bigSNP file: ", gt_get_file_names(x)[1])
    message("with backing file: ", gt_get_file_names(x)[2])
    message("make sure that you do NOT delete those files!")
    message("to reload the gen_tibble in another session, use:")
    message("gt_load('", file_name, "')")
  }
  return(c(file_name, gt_get_file_names(x)))
}


sub_gt <- function(path, replacement = "", stop_if_not_ext = TRUE) {
  pattern <- "\\.gt$"
  if (!grepl(pattern, path)) {
    stop("Path '%s' must have 'gt' extension.", path)
  }
  if (
    stop_if_not_ext &&
      (nchar(replacement) > 0) &&
      (substr(
        replacement,
        1,
        1
      ) !=
        ".")
  ) {
    stop("Replacement must be an extension starting with '.' if provided.")
  }
  sub(pattern, replacement, path)
}

#' Get the names of files storing the genotypes of a `gen_tibble`
#'
#' A function to return the names of the files used to store data in a
#' `gen_tibble`. Specifically, this returns the .rds file storing the big
#' @param x a [`gen_tibble`]
#' @returns a character vector with the names and paths of the two files
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # To retrieve the names of and paths to the .bk and .rds files use:
#' gt_get_file_names(example_gt)
#'
gt_get_file_names <- function(x) {
  if (inherits(x, "gen_tbl")) {
    x <- x$genotypes
  } else if (!inherits(x, "vctrs_bigSNP")) {
    stop("x should be a vctrs_bigSNP object")
  }

  return(c(
    # nolint
    attr(x, "bigsnp")$genotypes$rds,
    attr(x, "bigsnp")$genotypes$backingfile
  ))
}
