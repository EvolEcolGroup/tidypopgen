#' simple function to write a test plink BED file
#'
#' We start with data in the format that we expect to come out later on
#' @keywords internal

make_test_bed <- function(genotypes, loci, ind_meta, path_out){
  if (nrow(genotypes)!=nrow(loci)){
    genotypes <- t(genotypes)
  }

  bim <- genio::make_bim( n = nrow(genotypes) )
  bim$chr <- loci$chromosome
  bim$id <- loci$name
  bim$pos <- loci$position
  bim$alt <- loci$allele_alt
  bim$ref <- loci$allele_ref

  fam <- genio::make_fam(n= nrow(ind_meta))
  fam$fam <- ind_meta$population
  fam$id <- ind_meta$id
  rownames(genotypes) <- bim$id
  colnames(genotypes) <- fam$id

  genio::write_plink(path_out, genotypes, bim, fam)
}
