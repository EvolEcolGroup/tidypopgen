#' @rdname vcfR_conversion
#' @aliases vcfR2genlight
#'
#' @param n.cores integer specifying the number of cores to use.
#'
#' @examples
#' adegenet_installed <- require("adegenet")
#' if (adegenet_installed) {
#'   data(vcfR_test)
#'   # convert to genlight (preferred method with bi-allelic SNPs)
#'   gl <- vcfR2genlight(vcfR_test)
#'
#'   # convert to genind, keeping information about allelic state
#'   # (slightly slower, but preferred method for use with the "poppr" package)
#'   gid <- vcfR2genind(vcfR_test, return.alleles = TRUE)
#'
#'   # convert to genind, returning allelic states as 0, 1, 2, etc.
#'   # (not preferred, but slightly faster)
#'   gid2 <- vcfR2genind(vcfR_test, return.alleles = FALSE)
#' }
#'
#' @export
vcfR2gen_tibble <- function(x, n_cores=1){

  # subset to biallelic loci only
  bi <- is.biallelic(x)
  if(sum(!bi) > 0){
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), 'loci will be omitted from the genlight object.'))
    warning(msg)
    x <- x[bi,]
  }

  # fill in any missing IDs
  x <- vcfR::addID(x)

  # create loci table
  loci <- tibble(name = vcfR::getID(x),
                 chromosome = vcfR::getCHROM(x),
                 position = vcfR::getPOS(x),
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

  indiv_meta <- tibble(id = colnames(x), population = NA)

  x <- gen_tibble(indiv_meta= indiv_meta, genotypes = x, loci = loci)

  return(x)
}
