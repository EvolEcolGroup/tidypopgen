#' Constructor for a `gen_tibble`
#'
#' A `gen_tibble` stores genotypes for individuals in a tidy format. DESCRIBE
#' here the format
#' @param bigsnp_path the path to a [`bigsnpr::bigSNP`] file created with
#' [bigsnpr::snp_readBed()]
#' @returns an object of the class `gen_tbl`.
#' @export

gen_tibble <- function(bigsnp_path){
  if (is.character(bigsnp_path)){
    bigsnp_obj <- bigsnpr::snp_attach(bigsnp_path)
  } else {
    stop("bigsnp_path should be pointing to a bigsnp rds file")
  }
  ind_meta <- list(id = bigsnp_obj$fam$sample.ID,
                             population = bigsnp_obj$fam$family.ID)

  ind_meta$genotypes <- new_vctrs_bigsnp(bigsnp_obj, bigsnp_obj$fam$sample.ID)

  tibble::new_tibble(
    ind_meta,
    class = "gen_tbl"
  )
}



new_vctrs_bigsnp <- function(bigsnp_obj, names) {
  loci <- tibble::tibble(big_index = seq_len(nrow(bigsnp_obj$map)),
                         name = bigsnp_obj$map$marker.ID,
                         chromosome = bigsnp_obj$map$chromosome,
                         position = bigsnp_obj$map$physical.pos,
                         genetic_dist = bigsnp_obj$map$genetic.dist,
                         allele_ref = bigsnp_obj$map$allele2,
                         allele_alt = bigsnp_obj$map$allele1
  )
  vctrs::new_vctr(seq_len(nrow(bigsnp_obj$fam)),
                  bigsnp = bigsnp_obj,
                  loci=loci,
                  names=names,
                  class = "vctrs_bigSNP")
}

#' @export
summary.vctrs_bigSNP <- function(object, ...){
  summary(rep("bigSNP-genotypes",length(object)))
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

# print method
#' @export
tbl_sum.gen_tbl <- function(x, ...) {
  c(
    "A gen_tibble" = paste(nrow(show_loci(x))," loci"),
    NextMethod()
  )
}

