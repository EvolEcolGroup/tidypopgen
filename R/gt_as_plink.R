#' Export a `gen_tibble` object to PLINK bed format
#'
#' This function exports all the information of a `gen_tibble` object into
#' a PLINK bed file.
#'
#' @param x a [`gen_tibble`] object
#' @param bedfile a character string giving the path to output bed file.
#' @param overwrite boolean whether to overwrite the file.
#' @returns TRUE if successful
#' @export


gt_as_plink <- function(x, bedfile = NULL,
                           overwrite = TRUE){
  if (is.null(bedfile)){
    bedfile <- bigstatsr::sub_bk(attr(x$genotypes,"bigsnp")$genotypes$backingfile,".bed")
  }
  if (file_ext(bedfile)!="bed"){
    bedfile <- paste0(bedfile,".bed")
  }
  if (file.exists(bedfile)){
    if (overwrite){
      file.remove(bedfile)
      file.remove(bigsnpr::sub_bed(bedfile,".bim"))
      file.remove(bigsnpr::sub_bed(bedfile,".fam"))
    } else {
      stop(bedfile," already exists; remove if first or set 'overwrite' = TRUE")
    }
  }

  bed_path <- bigsnpr::snp_writeBed(attr(x$genotypes,"bigsnp"),
                        bedfile = bedfile,
                        ind.row = vctrs::vec_data(x$genotypes),
                        ind.col = show_loci(x)$big_index)
  # the bim and fam file only contain the original information in the bigSNP object
  # TODO we should update them with the info from the gentibble
  bed_path
}