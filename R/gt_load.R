#' Load a gen_tibble
#'
#' Load a `gen_tibble` previously saved with [gt_save()]. If the *.rds* and
#'  *.bk* files have not
#' been moved, they should be found automatically. If they were moved, use
#' `reattach_to` to point to the *.rds* file (the *.bk* file needs to be in the
#' same directory as the *.rds* file).
#'
#' @param file the file name, including the full path. If it does not end with
#'   *.gt*, the extension will be added.
#' @param reattach_to the file name, including the full path, of the *.rds* file
#'   if it was moved. It assumes that the *.bk* file is found in the same path.
#'   You should be able to leave this to NULL unless you have moved the files.
#' @returns a [`gen_tibble`]
#' @seealso [gt_save()]
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # remove some individuals
#' example_gt_filtered <- example_gt %>% filter(id != "a")
#'
#' # save the filtered gen_tibble object
#' backing_files <- gt_save(example_gt_filtered,
#'   file_name = paste0(tempfile(), "_example_filtered")
#' )
#'
#' # backing_files[1] contains the name of the saved .gt file
#' backing_files[1]
#'
#' # To load the saved gen_tibble object, use the path to the saved .gt file
#' reloaded_gt <- gt_load(backing_files[1])
#'
#' # And we have loaded the gt without individual "a"
#' reloaded_gt
gt_load <- function(file = NULL, reattach_to = NULL) {
  if (file_ext(file) != "gt") {
    file <- paste0(file, ".gt")
  }
  if (!file.exists(file)) {
    stop(file, " does not exist")
  }
  x <- readRDS(file)
  # update if we have old objects
  x <- update_old_bigsnp(x)

  # now deal with the backing files
  if (is.null(reattach_to)) {
    # we expect the backing files to exist
    backingfiles <- gt_get_file_names(x)
    if (!file.exists(backingfiles[1])) {
      stop(backingfiles[1], " does not exist")
    }
    if (!file.exists(backingfiles[2])) {
      stop(backingfiles[2], " does not exist")
    }
  } else {
    # we want to reattach to a backing file that has been moved
    if (file_ext(reattach_to) != "rds") {
      reattach_to <- paste0(reattach_to, ".rds")
    }
    # attach the new files
    attr(x$genotypes, "fbm") <- bigstatsr::big_attach(reattach_to)
    attr(x$genotypes, "fbm_md5sum") <- tools::md5sum(reattach_to)
  }

  return(x)
}


sub_rds <- function(path, replacement = "", stop_if_not_ext = TRUE) {
  pattern <- "\\.rds$"
  if (!grepl(pattern, path)) {
    stop(sprintf("Path '%s' must have '.rds' extension.", path))
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


# Update old objects to use FBM directly instead of a bigsnp_obj
update_old_bigsnp <- function(x) {
  # if we have a bignsp attribute, we need to upgrade
  if (!is.null(attr(x$genotypes, "bigsnp"))) {
    # move the fbm into its own slot
    attr(x$genotypes, "fbm") <- attr(x$genotypes, "bigsnp")$genotypes
    # fix ploidy
    if ((attr(x$genotypes, "ploidy") != 0) && (attr(x$genotypes, "ploidy") != -2)) { # nolint
      fbm_ploidy <- rep(attr(x$genotypes, "ploidy"), nrow(x))
    } else {
      fbm_ploidy <- attr(x, "bigsnp")$fam$ploidy[vctrs::vec_data(x)]
    }
    attr(x$genotypes, "fbm_ploidy") <- fbm_ploidy

    # save it (we have to overwrite the old file, as rds and bk have to have
    # the same name)
    saveRDS(attr(x$genotypes, "fbm"), file = gt_get_file_names(x)[1])
    # add the file to the object
    attr(x$genotypes, "fbm_file") <- gt_get_file_names(x)[1]
    # add the md5 sum of this file
    attr(x$genotypes, "fbm_md5sum") <- tools::md5sum(gt_get_file_names(x)[1])
    # now remove the old bignsp objects
    attr(x$genotypes, "bigsnp") <- NULL
    attr(x$genotypes, "bigsnp_file") <- NULL
    attr(x$genotypes, "bigsnp_md5sum") <- NULL
    message(
      "your gen_tibble was in an old format, and it has now ",
      "been updated\n",
      "note that the old RDS was overwritten"
    )
  }
  return(x)
}
