#' Save a gen_tibble
#'
#' Save the tibble (and update the backing files). The `gen_tibble` object is
#' saved to a file with extension *.gt*, together with update its *.rds* and
#' *.bk* files. Note that multiple *.gt* files can be linked to the same *.rds*
#' and *.bk* files; generally, this occurs when we create multiple subsets of
#' the data. The *.gt* file then stores the information on what subset of the
#' full dataset we are interested in, whilst the *.rds* and *.bk* file store the
#' full dataset. To reload a `gen_tibble`, you can pass the name of the *.gt*
#' file with [gt_load()]. The *.gt* file is compressed with gzip at its
#' fastest level, as saving the information on millions of loci with the
#' default level of [saveRDS()] is slow.
#'
#' If the names of the loci are integers (`names_as_int = TRUE`), their
#' original names are in a sidecar file next to the *.bk* file (*.loci.tsv*,
#' see [gt_lookup_loci()]); if the chromosomes are integers
#' (`chromosomes_as_int = TRUE`), their names are in another sidecar file
#' (*.chromosomes.tsv*, see [gt_lookup_chromosomes()]). The sidecars are
#' listed in the messages and in the value returned, after the *.bk* file.
#' Keep them with the other files (if you move the files, move them too, so
#' that [gt_load()] with `reattach_to` finds them).
#'
#' @param x a [`gen_tibble`]
#' @param file_name the file name, including the full path. If it does not end
#'   with *.gt*, the extension will be added.
#' @param quiet boolean to suppress information about the files
#' @returns the file name and path of the *.gt* file, followed by those of the
#'   *.rds* and *.bk* files and of any sidecar files (as given by
#'   [gt_get_file_names()])
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

  # we update the fbm object
  .gt_get_fbm(x)$save()

  if (is.null(file_name)) {
    file_name <- bigstatsr::sub_bk(gt_get_file_names(x)[2], ".gt")
  }
  if (file_ext(file_name) != "gt") {
    file_name <- paste0(file_name, ".gt")
  }
  # and now save our gen_tibble
  save_gt_file(x, file_name)
  if (!quiet) {
    message("\ngen_tibble saved to ", file_name)
    message("using FBM RDS: ", gt_get_file_names(x)[1])
    message("with FBM backing file: ", gt_get_file_names(x)[2])
    message_sidecar_files(x)
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
  save_gt_file(x, file_name)
  if (!quiet) {
    message("\ngen_tibble saved to ", file_name)
    message("using FBM RDS: ", gt_get_file_names(x)[1])
    message("with FBM backing file: ", gt_get_file_names(x)[2])
    message_sidecar_files(x)
    message("make sure that you do NOT delete those files!")
    message("to reload the gen_tibble in another session, use:")
    message("gt_load('", file_name, "')")
  }
  return(c(file_name, gt_get_file_names(x)))
}

#' Get the names of files storing the genotypes of a `gen_tibble`
#'
#' A function to return the names of the files used to store data in a
#' `gen_tibble`. Specifically, this returns the *.rds* file storing the big
#' matrix (FBM) and its *.bk* backing file, always as the first two elements
#' and in this order. If the `gen_tibble` has integer names of loci
#' (`names_as_int = TRUE`) or integer chromosomes (`chromosomes_as_int =
#' TRUE`), the sidecar files holding their original names (*.loci.tsv*, then
#' *.chromosomes.tsv*, see [gt_lookup_loci()] and [gt_lookup_chromosomes()])
#' follow. A sidecar is listed only if the `gen_tibble` uses it (its loci table
#' holds integer names, or integer chromosomes) and the file exists, so a
#' `gen_tibble` without integer names or chromosomes has exactly two files.
#' @param x a [`gen_tibble`]
#' @returns a character vector with the names and paths of the *.rds* and
#'   *.bk* files, followed by those of any sidecar files
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # To retrieve the names of and paths to the .bk and .rds files use:
#' gt_get_file_names(example_gt)
#'
gt_get_file_names <- function(x) {
  # the .rds and .bk files first (code indexes them with [1] and [2]), then
  # the sidecars that the gen_tibble uses and that exist (sidecar_files()
  # builds their paths with fbm_file_names(), not with this function)
  c(fbm_file_names(x), unname(sidecar_files(x)))
}

# the .rds and .bk files of the FBM of a gen_tibble (or of its genotypes
# column), without the sidecars
fbm_file_names <- function(x) {
  if (inherits(x, "gen_tbl")) {
    x <- x$genotypes
  } else if (!inherits(x, "vctrs_bigSNP")) {
    stop("x should be a vctrs_bigSNP object")
  } else if (is.null(attr(x, "fbm"))) {
    stop("x must be a genotype vector with attr(., 'fbm') present")
  }

  return(c(
    # nolint
    attr(x, "fbm")$rds,
    attr(x, "fbm")$backingfile
  ))
}

# save a gen_tibble to its .gt file with saveRDS(), compressed with gzip at
# level 1: compared with the default level of saveRDS() (6), this halves the
# time to save millions of loci, for a file about 20% larger, and the file
# is read by readRDS() as any .rds file
save_gt_file <- function(x, file_name) {
  file_con <- gzfile(file_name, open = "wb", compression = 1)
  on.exit(close(file_con))
  saveRDS(x, file_con)
}
