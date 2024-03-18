#' Export a `gen_tibble` object to PLINK ped format
#'
#' This function exports all the information of a `gen_tibble` object into
#' a PLINK ped file.
#'
#' This function generates a PLINK raw file with alternate allele information
#' (as outputted using `--include-alt` and
#' `--recodeA`), as well as a map file (as generated with `--recode`).
#'
#' @param x a [`gen_tibble`] object
#' @param file a character string giving the path to the file to convert,
#' with the extension ".raw" or "ped"
#' @param overwrite boolean whether to overwrite the file
#' @returns TRUE if successful
#' @export


gt_write_plink <- function(x, file = NULL,
                           overwrite = TRUE){
  if (is.null(file)){
    file <- bigstatsr::sub_bk(attr(x$genotypes,"bigsnp")$genotypes$backingfile,".bed")
  }
  if (file.exists(file)){
    if (overwrite){
      file.remove(file)
    } else {
      stop(file," already exists; remove if first or set 'overwrite' = TRUE")
    }
  }

  bed_path <- bigsnpr::snp_writeBed(attr(x$genotypes,"bigsnp"),
                        bedfile = file,
                        ind.row = vctrs::vec_data(x$genotypes),
                        ind.col = show_loci(x)$big_index)
  # the bim and fam file only contain the original information in the bigSNP object
  # we update them with the info from the gentibble
  bed_path
}
