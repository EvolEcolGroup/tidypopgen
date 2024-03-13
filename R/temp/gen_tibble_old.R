#' Constructor for a `gen_tibble`
#'
#' A `gen_tibble` stores genotypes for individuals in a tidy format. DESCRIBE
#' here the format
#'
#' @param ind_meta a list, data.frame or tibble with compulsory columns 'id'
#'  and 'population', plus any additional metadata of interest.
#' @param genotypes a matrix of counts of alternative alleles, one row per
#' individual and a column per locus
#' @param loci a data.frame or tibble, with compulsory columns 'name', 'chromosome',
#' and 'position'
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

gen_tibble <- function(genotypes, ind_meta = NULL, loci = NULL){
  # TODO check object types

  if (is.character(genotypes)){
    genotypes_big <- bigsnpr::snp_attach(genotypes)
  } else {
    stop("genotypes should be pointing to a bigsnpr rds file")
  }

  if (is.null(ind_meta)){
    ind_meta <- tibble::tibble(id = genotypes_big$fam$sample.ID,
                               population = genotypes_big$fam$family.ID)
  }
  if (!all(c("id", "population") %in% names(ind_meta))){
    stop("ind_meta does not include the compulsory columns 'id' and 'population")
  }

#  check_valid_loci(loci)
  ind_meta <- as.list(ind_meta)
  #TODO this could be parallelised
  ind_meta$genotypes <- vctrs_bigsnp(ind_meta$id, genotypes_big)

  browser()

  tibble::new_tibble(
    ind_meta,
    class = "gen_tbl"
  )
}



vctrs_bigsnp <- function(x = character(), bigsnp) {
  loci <- tibble::tibble(name = bigsnp$map$marker.ID,
                         chromosome = bigsnp$map$chromosome,
                         position = bigsnp$map$physical.pos,
                         genetic_dist = bigsnp$map$genetic.dist,
                         allele_ref = bigsnp$map$allele2,
                         allele_alt = bigsnp$map$allele1
  )
  vctrs::new_vctr(x, bigsnp = bigsnp, loci = loci, class = "vctrs_bigSNP")
}





check_valid_loci <- function(loci_df){
  loci_df <- as_tibble(loci_df)
  if (!all(c('name', 'chromosome', 'position','allele_ref','allele_alt') %in% names(loci_df))){
    stop("loci does not include the compulsory columns 'name', 'chromosome', 'position','allele_ref','allele_alt'")
  }
}



#' Test if a tibble is really `gen_tibble`
#'
#' Some `dplyr` operations strip the subclass from the tibble. This function
#' is used to check if the tibble is, in reality, still of class `gen_tbl`
#' @param .x the tibble
#' @returns NULL
#' @keywords internal

stopifnot_gen_tibble <- function(.x){
  if ("gentoypes" %in% names(.x)){
    stopifnot(.x$genotypes)
  }
}

#' # print method
#' #' @export
#' tbl_sum.gen_tbl <- function(x, ...) {
#'   c(
#'     "A gen_tibble" = paste(nrow(show_loci(x))," loci"),
#'     NextMethod()
#'   )
#' }

