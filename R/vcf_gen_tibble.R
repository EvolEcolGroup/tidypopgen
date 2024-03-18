#' Import VCF data into tidypopgen.
#'
#' Return  [`gen_tibble`] object.
#'
#' @param file a character path to the VCF file.
#' @param quiet whether to be verbose
#' @param ... further arguments to pass to [vcfR::read.vcfR()].
#' @returns a tibble with 3 elements: maf, missingness and hwe_p
#' @export

vcf_gen_tibble <- function(file, quiet = TRUE, ...) {
  x <- vcfR::read.vcfR(file = path, verbose = !quiet, ...)

  x <- vcfR::addID(x)

  # create loci table
  loci <- tibble(name = vcfR::getID(x),
                 chromosome = vcfR::getCHROM(x),
                 position = vcfR::getPOS(x),
                 #TODO: why is genetic_dist necessary here?
                 genetic_dist = 0,
                 allele_ref = vcfR::getREF(x),
                 allele_alt = vcfR::getALT(x))

  x <- vcfR::extract.gt(x)
  x[x=="0|0"] <- 0
  x[x=="0|1"] <- 1
  x[x=="1|0"] <- 1
  x[x=="1|1"] <- 2
  x[x=="0/0"] <- 0
  x[x=="0/1"] <- 1
  x[x=="1/0"] <- 1
  x[x=="1/1"] <- 2
  # make sure these are numeric
  x <- apply(x, 2, as.numeric)

  ind_meta <- tibble(id = colnames(x), population = NA)

  # using the gen_tibble.matrix method
  gen_tibble(x = x, indiv_meta = ind_meta, loci = loci)
}
