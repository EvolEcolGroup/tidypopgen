#' Simple function to write a plink BED file from data.frames
#'
#' This function write a bed file using information in tables with the same format
#' as they are stored in various parts of a [`gen_tibble`]. This is mostly used
#' for creating small test files, but could also be used to convert small datasets
#' manually into BED files (and the appropriate BIM and FAM companions).
#'
#' @param genotypes a matrix of genotype dosages
#' @param loci a table of loci information (as used in [`gen_tibble`])
#' @param indiv_meta a table with 'id' and 'population'
#' @param path_out is the path and the prefix, to which '.bed', '.bim', and '.fam'
#' will be added.
#' @returns NULL it writes the files
#' @export

gt_write_bed_from_dfs <- function(genotypes, loci, indiv_meta, path_out=NULL){

  if (is.null(path_out)){
    stop ("path_out needs to be given")
  }
  if (nrow(genotypes)!=nrow(loci)){
    genotypes <- t(genotypes)
  }

  bim <- genio::make_bim( n = nrow(genotypes) )
  bim$chr <- loci$chromosome
  bim$id <- loci$name
  bim$pos <- loci$position
  # bim$posg <- loci$genetic_dist
  bim$ref <- loci$allele_ref
  bim$alt <- loci$allele_alt
  # replace missing alleles with 0 (the plink default)
  bim$alt[is.na(bim$alt)] <- "0"

  fam <- genio::make_fam(n= nrow(indiv_meta))
  fam$fam <- indiv_meta$population
  fam$id <- indiv_meta$id
  rownames(genotypes) <- bim$id
  colnames(genotypes) <- fam$id

  genio::write_plink(path_out, genotypes, bim, fam, verbose=FALSE)
  return(paste0(path_out,".bed"))
}
