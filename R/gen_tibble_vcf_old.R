# read in a vcf
gen_tibble_vcf2 <- function(x, ...,
                           valid_alleles = c("A", "T", "C", "G"),
                           missing_alleles = c("0","."),
                           backingfile = NULL, quiet = FALSE) {
  x <- vcfR::read.vcfR(file = x, verbose = !quiet, ...)

  x <- vcfR::addID(x)

  # create loci table
  loci <- tibble(name = vcfR::getID(x),
                 chromosome = vcfR::getCHROM(x),
                 position = vcfR::getPOS(x),
                 genetic_dist = 0,
                 allele_ref = vcfR::getREF(x),
                 allele_alt = vcfR::getALT(x))

  bi <- vcfR::is.biallelic(x)
  x <- vcfR::extract.gt(x)
  x <- x[bi,]
  # TODO from the first locus, we should figure out the ploidy
  # for the moment, we hardcode to ploidy 2
  ploidy = 2

  x[x=="0|0"] <- 0
  x[x=="0|1"] <- 1
  x[x=="1|0"] <- 1
  x[x=="1|1"] <- 2
  x[x=="0/0"] <- 0
  x[x=="0/1"] <- 1
  x[x=="1/0"] <- 1
  x[x=="1/1"] <- 2
  # additional conversion for tetraploids
  # but it would be better to make the conversion more generic for all ploidies
  # x[x == "1/1/1/1"] <- 4
  # x[x == "0/1/1/1"] <- 3
  # x[x == "0/0/1/1"] <- 2
  # x[x == "0/0/0/1"] <- 1
  # x[x == "0/0/0/0"] <- 0
  # make sure these are numeric
  x <- apply(x, 2, as.numeric)

  ind_meta <- tibble(id = colnames(x), population = NA)

  # using the gen_tibble.matrix method
  new_gen_tbl <- gen_tibble(x = t(x),
                            indiv_meta = ind_meta,
                            loci = loci,
                            backingfile = backingfile,
                            ploidy = ploidy)
  check_allele_alphabet (new_gen_tbl, valid_alleles = valid_alleles,
                         missing_alleles = missing_alleles)
  show_loci(new_gen_tbl) <- harmonise_missing_values(show_loci(new_gen_tbl), missing_alleles = missing_alleles)
  return(new_gen_tbl)

}
