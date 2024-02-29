#' Constructor for a `gen_tibble`
#'
#' A `gen_tibble` stores genotypes for individuals in a tidy format. DESCRIBE
#' here the format
#' @param ind_meta a list, data.frame or tibble with compulsory columns 'id'
#'  and 'population', plus any additional metadata of interest.
#' @param genotypes a matrix of counts of alternative alleles, one row per
#' individual and a column per locus
#' @param loci a data.frame or tibble, with compulsory columns 'name', 'chromosome',
#' and 'position'
#' @param ploidy a vector giving the ploidy of each individual. If 'ploidy' is not
#' specified, diploids are assumed.
#' @returns an object of the class `gen_tbl`.
#' @examples
#' test_ind_meta <- data.frame (id=c("a","b","c"),
#'                              population = c("pop1","pop1","pop2"))
#' test_genotypes <- rbind(c(1,1,0,1,1,0),
#'                         c(2,1,0,0,0,0),
#'                         c(2,2,0,0,1,1))
#' test_loci <- data.frame(name=paste0("rs",1:6),
#'                         chromosome=c(1,1,1,1,2,2),
#'                         position=c(3,5,65,343,23,456),
#'                         allele_ref = c("a","t","c","g","c","t"),
#'                         allele_alt = c("t","c", NA,"c","g","a"))
#' test_gen <- gen_tibble(ind_meta = test_ind_meta,
#'                        genotypes = test_genotypes,
#'                        loci = test_loci)
#' test_gen
#' @export

gen_tibble <- function(ind_meta, genotypes, loci, ploidy=NULL){
  # TODO check object types
  if (!all(c("id", "population") %in% names(ind_meta))){
    stop("ind_meta does not include the compulsory columns 'id' and 'population")
  }
  if (is.null(ploidy)){
    ploidy <- rep(2, nrow(ind_meta))
  }
  if (!all(c('name', 'chromosome', 'position','allele_ref','allele_alt') %in% names(loci))){
    stop("loci does not include the compulsory columns 'name', 'chromosome', 'position','allele_ref','allele_alt'")
  }
  ind_meta <- as.list(ind_meta)
  #browser()
  ind_meta$genotypes <- lapply(1:nrow(genotypes), function(i) methods::new("SNPbin",
                                                                           as.integer(genotypes[i,]),ploidy=ploidy[i]) )
  attr(ind_meta$genotypes,"loci")<-tibble::as_tibble(loci)

  tibble::new_tibble(
    ind_meta,
    class = "gen_tbl"
  )
}
