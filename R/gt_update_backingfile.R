#' Update the backing matrix
#'
#' This functions forces a re-write of the file backing matrix to match the
#' [`gen_tibble`]. Individuals and loci are subsetted and reordered according to
#' the current state of the `gen_tibble`. A `.gt`, `.bk` and `.rds` file will be
#' created.
#'
#' This function does not check whether the positions of your genetic loci are
#' sorted. To check this, and update the file backing matrix, use
#' `gt_order_loci()`. Tests for this function are in
#' test_gt_order_loci.R
#'
#' @param .x a `gen_tibble` object
#' @param backingfile the path, including the file name without extension, for
#'   backing files used to store the data (they will be given a .bk and .rds
#'   automatically). If left to NULL (the default), the file name will be based
#'   on the name of the current backing file.
#' @param rm_unsorted_dist boolean to set `genetic_dist` to zero (i.e. remove
#'   it) if it is unsorted within the chromosomes.
#' @param chunk_size the number of loci to process at once
#' @param quiet boolean to suppress information about the files
#' @returns a [`gen_tibble`] with updated `.gt`, `.bk`, and `.rds` files
#' (i.e. a new File Backed Matrix)
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Here, new backingfiles are created, but without using `<-` the object
#' # loaded in the R session is not updated
#' gt_update_backingfile(example_gt)
#'
#' # Make sure to use `<-` to update the file names associated with the
#' # gen_tibble object loaded in the R session
#' example_gt <- example_gt %>% gt_update_backingfile()
#'
# Note that this is tested through gt_order_loci()
gt_update_backingfile <- function(
  .x,
  backingfile = NULL,
  chunk_size = NULL,
  rm_unsorted_dist = TRUE,
  quiet = FALSE
) {
  # if the backingfile is null, create a name based on the current backing file
  if (is.null(backingfile)) {
    backingfile <- change_duplicated_file_name(gt_get_file_names(.x)[2])
  }
  # if the chunk size is null, set it to the block size
  if (is.null(chunk_size)) {
    chunk_size <- bigstatsr::block_size(nrow(.x))
  }

  # set up chunks
  no_variants <- count_loci(.x)
  chunks_vec <- c(
    rep(chunk_size, floor(no_variants / chunk_size)),
    no_variants %% chunk_size
  )
  # remove any 0 length chunks (there could be one at the end)
  chunks_vec <- chunks_vec[chunks_vec > 0L]
  chunks_vec <- c(0, cumsum(chunks_vec))

  # get the 256 code
  original_256code <- attr(.x$genotypes, "fbm")$code256
  original_imputation_code <- attr(.x$genotypes, "imputed")

  # initialise a FBM with the dimensions of the data used in the gen_tibble
  new_bk_matrix <- bigstatsr::FBM.code256(
    nrow = length(.gt_fbm_rows(.x)),
    ncol = length(.gt_fbm_cols(.x)),
    code = original_256code,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )

  # now copy the raw values of the genotypes (so that imputed values are kept)
  # in chunks of loci; the copy is done in C++, so the genotypes are never
  # decoded into R
  fbm_rows <- .gt_fbm_rows(.x)
  fbm_cols <- .gt_fbm_cols(.x)
  for (i in seq_along(chunks_vec[-1])) {
    this_chunk <- seq(chunks_vec[i] + 1, chunks_vec[i + 1])
    gt_copy_fbm(
      BM_from = .gt_get_fbm(.x),
      BM_to = new_bk_matrix,
      rowInd = fbm_rows,
      colInd = fbm_cols[this_chunk],
      col_offset = chunks_vec[i],
      row_offset = 0,
      swap = logical(0),
      swap_lut = 0:255
    )
  }
  # get fbm_ploidy
  fbm_ploidy <- attr(.x$genotypes, "fbm_ploidy")
  # if it is not null, then subset it to the new individuals
  if (!is.null(fbm_ploidy)) {
    fbm_ploidy <- fbm_ploidy[.gt_fbm_rows(.x)]
  }

  # if we remove unsorted genetic distance, set it to zero if it is not sorted
  if (rm_unsorted_dist) {
    # check the genetic distances within each chromosome (splitting only that
    # column, rather than the whole loci table)
    if (any_unsorted_by_group(
      show_loci(.x)$genetic_dist,
      cast_chromosome_to_int(show_loci(.x)$chromosome)
    )) {
      if (!quiet) {
        message("Genetic distances are not sorted, setting them to zero")
      }
      show_loci(.x)$genetic_dist <- 0
    }
  }

  # save it
  fbm_file <- bigstatsr::sub_bk(new_bk_matrix$backingfile, ".rds")
  saveRDS(new_bk_matrix, fbm_file)

  # now return a gen_tibble based on the new backing file
  new_gen_tbl <- .x

  new_gen_tbl$genotypes <- new_vctrs_bigsnp(
    fbm_obj = new_bk_matrix,
    fbm_file = fbm_file,
    loci = show_loci(.x),
    indiv_id = .x$id,
    ploidy = attr(.x$genotypes, "ploidy"),
    fbm_ploidy = fbm_ploidy
  )

  # reorder the bignsp column in the loci table
  show_loci(new_gen_tbl)$big_index <- 1:count_loci(new_gen_tbl)

  # derive .gt filename from .bk backing file path
  file_name <- bigstatsr::sub_bk(gt_get_file_names(new_gen_tbl)[2], ".gt")

  # add imputation code
  if (!is.null(original_imputation_code)) {
    attr(new_gen_tbl$genotypes, "imputed") <- original_imputation_code
  }

  gt_save_light(new_gen_tbl, file_name = file_name, quiet = TRUE)

  if (!quiet) {
    message("\ngen_backing files updated, now")
    message("using FBM RDS: ", gt_get_file_names(new_gen_tbl)[1])
    message("with FBM backing file: ", gt_get_file_names(new_gen_tbl)[2])
    message("make sure that you do NOT delete those files!")
    message("to reload the gen_tibble in another session, use:")
    message("gt_load('", file_name, "')")
  }

  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(new_gen_tbl)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(new_gen_tbl) <- obj_class
  }

  return(new_gen_tbl)
}
