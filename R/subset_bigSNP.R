#' Subset a bigSNP object
#'
#' This function subsets an [`bigsnpr::bigSNP`] object, creating a new backing
#' file with the information subsetted and reorders as required.
#' @param X the [`bigsnpr::bigSNP`] object to subset
#' @param indiv_indices the indices of the individuals (rows) to keep
#' @param loci_indices the indices of the loci (columns) to keep
#' @param swap_indices the indices of loci where the alleles should be swapped
#'   (both in the bim table and in the genotypes)
#' @param backingfile the backing file (if null, a tempfile will be used)
#' @returns a [`bigsnpr::bigSNP`] object
#' @keywords internal
#' @noRd
# nolint start
subset_bigSNP <- function(
    X, # nolint end
    indiv_indices = NULL,
    loci_indices = NULL,
    swap_indices = NULL,
    backingfile = NULL) {
  if (is.null(indiv_indices)) {
    indiv_indices <- bigstatsr::rows_along(X)
  }
  if (is.null(loci_indices)) {
    loci_indices <- bigstatsr::cols_along(X)
  }
  new_map <- X$map[loci_indices, ]
  new_map_ref <- new_map # a copy to use when swapping alleles

  # save a copy of the FBM with only the columns and rows needed, in order
  if (is.null(backingfile)) {
    backingfile <- tempfile()
  }
  new_geno <- bigstatsr::big_copy(
    X$genotypes,
    ind.row = indiv_indices,
    ind.col = loci_indices,
    backingfile = backingfile
  )

  # now we swap the loci (if necessary)
  if (!is.null(swap_indices)) {
    # convert swap indices to the new order
    swap_indices <- match(swap_indices, loci_indices)
    swap_indices <- swap_indices[!is.na(swap_indices)]

    # store the old code
    code_ref <- new_geno$code256
    # set a code that allows us to see both raw data and imputed data
    new_geno$code256 <- 1:256
    swap_locus <- function(x) {
      # note that we work on the raw bytes valus
      x_new <- x
      x_new[x == 3] <- 1
      x_new[x == 1] <- 3
      x_new[x == 5] <- 7
      x_new[x == 7] <- 5
      as.raw(x_new - 1)
    }
    # change the genotypes
    for (i in swap_indices) {
      new_geno[, i] <- swap_locus(new_geno[, i])
    }
    # now swap the alleles in the map
    new_map$allele1[swap_indices] <- new_map_ref$allele2[swap_indices]
    new_map$allele2[swap_indices] <- new_map_ref$allele1[swap_indices]
    # now reset the code
    new_geno$code256 <- code_ref
  }

  # now reassemble
  # Create the bigSNP object
  snp_list <- structure(
    list(
      genotypes = new_geno,
      fam = X$fam[indiv_indices, ],
      map = new_map
    ),
    class = "bigSNP"
  )

  # save it and return the path of the saved object
  rds <- bigstatsr::sub_bk(new_geno$backingfile, ".rds")
  saveRDS(snp_list, rds)
  # we return the object
  snp_list
}
