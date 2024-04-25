#' Read loci information directly from PLINK raw files
#'
#' This function extract locus names and alleles from a PLINK raw file, without
#' the need to build a full `gen_tibble` object. This function returns a tibble
#' equivalent to a subset of `loci` attribute of the `genotypes` column (i.e
#' missing 'chromosome' and 'position')
#' @param file the path to the PLINK raw file
#' @returns an tibble with the loci inofrmation
#' @keywords internal

read_plink_raw_loci <- function(file){
  # we use the same code as in adegenet::read.PLINK
  ext <- adegenet::.readExt(file)
  ext <- toupper(ext)
  if(ext != "RAW") warning("wrong file extension - '.raw' expected")
  col_names <- scan(file,what="character",sep=" ",quiet=TRUE,  nlines=1, blank.lines.skip=FALSE)
  indiv_meta_list <- lapply(1:6,function(i) NULL)
  names(indiv_meta_list) <- col_names[1:6]
  loci_names <- col_names[7:length(col_names)]
  # remove underscore followed by a digit at the end of locus name
  # why would that be a problem?!?
  #  loci_names <- gsub("_[1-9]$","",loci_names)

  # now parse the alleles
  #TODO we should check that we have info for both alleles

  # parse the information from the loci col names
  # the format is rs3131972_A(/G)
  loci <- tibble::tibble(name = substr(loci_names,1,nchar(loci_names)-6),
                         allele_ref = tolower(substr(loci_names,nchar(loci_names)-1,nchar(loci_names)-1)),
                         allele_alt = tolower(substr(loci_names,nchar(loci_names)-4,nchar(loci_names)-4)))
  # replace zeroes (missing) with NA
  loci$allele_alt[loci$allele_alt=="0"] <- NA
  return(loci)
}
