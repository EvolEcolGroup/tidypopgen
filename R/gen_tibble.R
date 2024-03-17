#' Constructor for a `gen_tibble`
#'
#' A `gen_tibble` stores genotypes for individuals in a tidy format. DESCRIBE
#' here the format
#' @param x can be:
#' - a string giving the path to a plink BED file, or a [`bigsnpr::bigSNP`]
#' - a string giving the path to a RDS file storing a `bigSNP` object from
#' the `bigsnpr` package (usually created with [bigsnpr::snp_readBed()])
#' - a genotype matrix of dosages (0, 1, 2, NA) giing the dosage of the alternate
#' allele.
#' @param ... unused (necessary for the method to accept different parameters depending
#' on `x`)
#' @param backingfile the path, including the file name without extension,
#' for backing files used to store the data (they will be given a .bk
#' and .RDS automatically). This is not needed if `x` is already an .RDS file.
#' If `x` is a .BED file and `backingfile` is left NULL, the backing file will
#' be saved in the same directory as the
#' bed file, using the same file name but with a different file type (.bk rather
#' than .bed). If `x` is a genotype matrix and `backingfile` is NULL, then a
#' temporary file will be created (but note that R will delete it at the end of
#' the session!)
#' @param quiet provide information on the files used to store the data
#' @returns an object of the class `gen_tbl`.
#' @rdname gen_tibble
#' @export

gen_tibble <- function(x, ...,  backingfile = NULL, quiet = FALSE) {
  UseMethod("gen_tibble", x)
}

#' @export
#' @rdname gen_tibble
gen_tibble.character <- function(x, ..., backingfile = NULL, quiet = FALSE){
  rlang::check_dots_empty()
  # if it is a bed file, we convert it to a bigsnpr
  if (tolower(file_ext(x))=="bed"){
    if (is.null(backingfile)){
      backingfile <- bigsnpr::sub_bed(x)
    }
    bigsnp_path <- bigsnpr::snp_readBed(x,
                                        backingfile = backingfile)
  }  else if (tolower(file_ext(x))=="rds"){
    bigsnp_path <- x
  } else {
    stop("file_path should be pointing to a either a PLINK .bed file or a bigSNP .rds file")
  }

  bigsnp_obj <- bigsnpr::snp_attach(bigsnp_path)
  if (!quiet){
    message("\n\nusing bigSNP file: ", bigsnp_path)
    message("with backing file: ", bigsnp_obj$genotypes$backingfile)
    message("make sure that you keep those files and don't delete them!")
  }

  indiv_meta <- list(id = bigsnp_obj$fam$sample.ID,
                             population = bigsnp_obj$fam$family.ID)

  indiv_meta$genotypes <- new_vctrs_bigsnp(bigsnp_obj,
                                           bigsnp_file = bigsnp_path,
                                           indiv_id = bigsnp_obj$fam$sample.ID)

  tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )
}


# simple function to extract the extension of a file
file_ext <- function(x){
  utils::tail(unlist(strsplit(x,".",fixed = TRUE)),n=1)
}

#' create a vctrs_bigSNP
#' @param bigsnp_obj the bigsnp object
#' @param bigsnp_file the file to which the bigsnp object was saved
#' @param indiv_id ids of individuals
#' @returns a vctrs_bigSNP object
#' @keywords internal
new_vctrs_bigsnp <- function(bigsnp_obj, bigsnp_file, indiv_id) {
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
                  bigsnp_file = bigsnp_file,
                  loci=loci,
                  names=indiv_id,
                  class = "vctrs_bigSNP")
}

#' @export
summary.vctrs_bigSNP <- function(object, ...){
  summary(rep("bigSNP-genotypes",length(object)))
}


check_valid_loci <- function(loci_df){
  loci_df <- as_tibble(loci_df)
  if (!all(c('name', 'chromosome', 'position','genetic_dist', 'allele_ref','allele_alt') %in% names(loci_df))){
    stop("loci does not include the compulsory columns 'name', 'chromosome', 'position','allele_ref','allele_alt'")
  }
}


#' Create a bigSNP object from data.frames
#'
#' This function expects the indiv_meta and loci to have the correct columns
#' @param genotypes a genotype matrix
#' @param indiv_meta the individual meta information
#' @param loci the loci table
#' @keywords internal
bignsnp_from_dfs <- function(genotypes, indiv_meta, loci){
  warning("this function still needs testing")
  check_valid_loci(loci)
  bigGeno[] <- genotypes
  fam <- tibble(family.ID = indiv_meta$population,
                sample.ID = indiv_meta$id,
                paternal.ID = 0,
                maternal.ID = 0,
                sex = 0,
                affection = 0)
  map <- tibble(chromosome = loci$chromosome,
                marker.ID = loci$name,
                genetic.dist = loci$genetic_dist,
                physical.pos = loci$position,
                allele1 = loci$allele_alt,
                allele2 = loci$allele_ref)
  # Create the bigSNP object
  snp.list <- structure(list(genotypes = bigGeno, fam = fam, map = map),
                        class = "bigSNP")

  # save it and return the path of the saved object
  rds <- bigstatsr::sub_bk(bigGeno$backingfile, ".rds")
  saveRDS(snp.list, rds)
}







#################################################################################

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


bignsnp_from_dfs <- function(genotypes, indiv_meta, loci){
  bigGeno[] <- genotypes
  fam <- tibble(family.ID = indiv_meta$population,
                sample.ID = indiv_meta$id,
                paternal.ID = 0,
                maternal.ID = 0,
                sex = 0,
                affection = 0)
  map <- tibble(chromosome = loci$chromosome,
                marker.ID = loci$name,
                genetic.dist = loci$genetic_dist,
                physical.pos = loci$position,
                allele1 = loci$allele_alt,
                allele2 = loci$allele_ref)
  # Create the bigSNP object
  snp.list <- structure(list(genotypes = bigGeno, fam = fam, map = bim),
                        class = "bigSNP")

  # save it and return the path of the saved object
  rds <- sub_bk(bigGeno$backingfile, ".rds")
  saveRDS(snp.list, rds)
}
