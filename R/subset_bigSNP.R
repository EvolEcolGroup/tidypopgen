#' Subset a bigSNP object
#'
#' This function subsets an [`bigsnpr::bigSNP`] object, creating a new backing file with the
#' information subsetted and reorders as required.
#' @param X the [`bigsnpr::bigSNP`] object to subset
#' @param indiv_indices the indices of the individuals (rows) to keep
#' @param loci_indices the indices of the loci (columns) to keep
#' @param swap_indices the indices of loci where the alleles should be swapped (both
#' in the bim table and in the genotypes)
#' @param backingfile the backing file (if null, a tempfile will be used)
#' @returns a [`bigsnpr::bigSNP`] object
#' @keywords internal

subset_bigSNP <- function(X, indiv_indices=NULL, loci_indices=NULL, swap_indices=NULL,
                          backingfile = NULL){
  if (is.null(indiv_indices)){
    indiv_indices <- bigstatsr::rows_along(X)
  }
  if (is.null(loci_indices)){
    loci_indices <- bigstatsr::cols_along(X)
  }
  new_map <- X$map
  # first we swap the loci (if necessary)
  if (!is.null(swap_indices)){
    swap_locus <- function(x){
      x_new <- x
      x_new[x==2]<-0
      x_new[x==0]<-2
      x_new[is.na(x_new)]<-3
      x_new
    }
    # change the genotypes
    for (i in swap_indices){
      X$genotypes[,i]<-swap_locus(X$genotypes[,i])
    }
    # now swap the alleles in the map
    new_map$allele1[swap_indices]<-X$map$allele2[swap_indices]
    new_map$allele2[swap_indices]<-X$map$allele1[swap_indices]
  }

  if (is.null(backingfile)){
    backingfile <- tempfile(tmpdir = getOption("FBM.dir"))
  }
  # save a copy of the FBM with only the columns and rows needed, in the order necessary
  new_geno <- bigstatsr::big_copy(X$genotypes, ind.row = indiv_indices,
                      ind.col = loci_indices, backingfile = backingfile)
  # now reassemble
  # Create the bigSNP object
  snp.list <- structure(list(genotypes = new_geno,
                             fam = X$fam[indiv_indices, ],
                             map = new_map[loci_indices, ]),
                            class = "bigSNP")

  # save it and return the path of the saved object
  rds <- bigstatsr::sub_bk(new_geno$backingfile, ".rds")
  saveRDS(snp.list, rds)
  # we return the object
  snp.list
}

