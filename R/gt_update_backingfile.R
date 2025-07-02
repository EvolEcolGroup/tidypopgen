#' Update the backing matrix
#'
#' This functions forces a re-write of the file backing matrix to match the
#' [`gen_tibble`]. Individuals and loci are subsetted and reordered according to
#' the current state of the `gen_tibble`. Tests for this function are in
#' test_gt_order_loci.R
#'
#' This function does not check whether the positions of your genetic loci are
#' sorted. To check this, and update the file backing matrix, use
#' `gt_order_loci()`.
#'
#' @param .x a `gen_tibble` object
#' @param backingfile the path, including the file name without extension, for
#'   backing files used to store the data (they will be given a .bk and .RDS
#'   automatically). If left to NULL (the default), the file name will be based
#'   on the name f the current backing file.
#' @param rm_unsorted_dist boolean to set `genetic_dist` to zero (i.e. remove
#'   it) if it is unsorted within the chromosomes.
#' @param chunk_size the number of loci to process at once
#' @param quiet boolean to suppress information about the files
#' @returns a [`gen_tibble`] with a backing file (i.e. a new File Backed Matrix)
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% gt_update_backingfile()
#'
gt_update_backingfile <- function(
    .x,
    backingfile = NULL,
    chunk_size = NULL,
    rm_unsorted_dist = TRUE,
    quiet = FALSE) {
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
  chunks_vec <- c(0, cumsum(chunks_vec))

  # get the 256 code
  original_256code <- attr(.x$genotypes, "bigsnp")$genotypes$code256
  on.exit(
    attr(.x$genotypes, "bigsnp")$genotypes$code256 <- original_256code
  )

  # initialise a FBM with the dimensions of the data used in the gen_tibble
  new_bk_matrix <- bigstatsr::FBM.code256(
    nrow = length(.gt_bigsnp_rows(.x)),
    ncol = length(.gt_bigsnp_cols(.x)),
    code = original_256code,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )

  # set the code of original file as 1:256 to get raw values
  attr(.x$genotypes, "bigsnp")$genotypes$code256 <- as.integer(0:255)
  # now process the genotype matrix in chunks
  for (i in seq_along(chunks_vec[-1])) {
    this_chunk <- seq(chunks_vec[i] + 1, chunks_vec[i + 1])
    this_genotypes <- show_genotypes(.x, loci_indices = this_chunk)
    new_bk_matrix[, this_chunk] <- as.raw(this_genotypes)
  }
  # get the fam and map tables from the bigsnpr object
  fam <- attr(.x$genotypes, "bigsnp")$fam
  # reorder according to the gen_tibble
  fam <- fam[.gt_bigsnp_rows(.x), ]
  # same for map
  map <- attr(.x$genotypes, "bigsnp")$map
  map <- map[.gt_bigsnp_cols(.x), ]
  # if we remove unsorted genetic distance, set it to zero if it is not sorted
  if (rm_unsorted_dist) {
    if (
      any(unlist(
        show_loci(.x) %>%
          group_by(.data$chr_int) %>% # nolint start
          group_map(~ is.unsorted(.x$genetic_dist))
      )) ||
        any(unlist(
          show_loci(.x) %>%
            group_by(.data$chr_int) %>%
            group_map(~ duplicated(.x$genetic_dist))
        ))
    ) {
      # nolint end
      if (!quiet) {
        message("Genetic distances are not sorted, setting them to zero")
      }
      show_loci(.x)$genetic_dist <- 0
      map$genetic.dist <- 0
    }
  }

  # Create the bigSNP object
  bigsnp_obj <- structure(
    list(genotypes = new_bk_matrix, fam = fam, map = map),
    class = "bigSNP"
  )

  # save it
  bigsnp_file <- bigstatsr::sub_bk(new_bk_matrix$backingfile, ".rds")
  saveRDS(bigsnp_obj, bigsnp_file)

  # now return a gen_tibble based on the new backing file
  new_gen_tbl <- .x

  # now create a genotype column with this data
  new_gen_tbl$genotypes <- vctrs::new_vctr(
    seq_len(nrow(bigsnp_obj$fam)),
    bigsnp = bigsnp_obj,
    # TODO is this redundant with the info in the bigSNP object?
    bigsnp_file = bigsnp_file,
    # TODO make sure this does not take too long
    bigsnp_md5sum = tools::md5sum(bigsnp_file),
    loci = show_loci(.x),
    names = .x$id,
    ploidy = attr(.x$genotypes, "ploidy"),
    class = "vctrs_bigSNP"
  )
  # reorder the bignsp column in the loci table
  show_loci(new_gen_tbl)$big_index <- 1:count_loci(new_gen_tbl)

  if (!quiet) {
    message("\ngen_backing files updated, now")
    message("using bigSNP file: ", gt_get_file_names(new_gen_tbl)[1])
    message("with backing file: ", gt_get_file_names(new_gen_tbl)[2])
    message("make sure that you do NOT delete those files!")
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
