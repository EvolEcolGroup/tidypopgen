#' Constructor for a `gen_tibble`
#'
#' A `gen_tibble` stores genotypes for individuals in a tidy format. DESCRIBE
#' here the format
#' @param x can be:
#' - a string giving the path to a plink BED file, or a [`bigsnpr::bigSNP`]
#' - a string giving the path to a RDS file storing a `bigSNP` object from
#' the `bigsnpr` package (usually created with [bigsnpr::snp_readBed()])
#' - a string giving the path to a vcf file. Note that we currently read the whole
#' vcf in memory with `vcfR`, so only smallish vcf can be imported. Only biallelic
#' SNPs will be considered.
#' - a genotype matrix of dosages (0, 1, 2, NA) giing the dosage of the alternate
#' allele.
#' @param indiv_meta a list, data.frame or tibble with compulsory columns 'id'
#'  and 'population', plus any additional metadata of interest.
#' @param loci a data.frame or tibble, with compulsory columns 'name', 'chromosome',
#' and 'position','genetic_dist', 'allele_ref' and 'allele_alt'
#' @param ... if `x` is the name of a vcf file, additional arguments
#' passed to [vcfR::read.vcfR()]. Otherwise, unused.
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

  if ((tolower(file_ext(x))=="bed") || (tolower(file_ext(x))=="rds")){
    rlang::check_dots_empty()
    gen_tibble_bed_rds(x = x, ..., backingfile = backingfile, quiet = quiet)
  } else if (tolower(file_ext(x))=="vcf"){
    gen_tibble_vcf(x = x, ..., backingfile = backingfile, quiet = quiet)
  } else {
    stop("file_path should be pointing to a either a PLINK .bed file or a bigSNP .rds file")
  }
}


gen_tibble_bed_rds <- function(x, ..., backingfile = NULL, quiet = FALSE){

  # if it is a bed file, we convert it to a bigsnpr
  if (tolower(file_ext(x))=="bed"){
    if (is.null(backingfile)){
      backingfile <- bigsnpr::sub_bed(x)
    }
    bigsnp_path <- bigsnpr::snp_readBed(x,
                                        backingfile = backingfile)
  }  else if (tolower(file_ext(x))=="rds"){
    bigsnp_path <- x
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

gen_tibble_vcf <- function(x, ..., backingfile = NULL, quiet = FALSE) {
  x <- vcfR::read.vcfR(file = x, verbose = !quiet, ...)

  x <- vcfR::addID(x)

  # create loci table
  loci <- tibble(name = vcfR::getID(x),
                 chromosome = vcfR::getCHROM(x),
                 position = vcfR::getPOS(x),
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
  gen_tibble(x = x,
             indiv_meta = ind_meta,
             loci = loci,
             backingfile = backingfile)
}

#' @export
#' @rdname gen_tibble
gen_tibble.matrix <- function(x, indiv_meta, loci, ..., backingfile = NULL, quiet = FALSE){
  rlang::check_dots_empty()
  # TODO check object types
  if (!all(c("id", "population") %in% names(indiv_meta))){
    stop("ind_meta does not include the compulsory columns 'id' and 'population")
  }
  # check that x (the genotypes) is numeric matrix
  if (inherits(x,"data.frame")){
    x <- as.matrix(x)
  }
  if (any(!inherits(x,"matrix"), !is.numeric(x))){
    stop("'x' is not a matrix of integers")
  }

  # use code for NA in FBM.256
  x[is.na(x)]<-3


  bigsnp_obj <- gt_write_bigsnp_from_dfs(genotypes = x,
                                          indiv_meta = indiv_meta,
                                          loci = loci,
                                          backingfile = backingfile)
  bigsnp_path <- bigstatsr::sub_bk(bigsnp_obj$genotypes$backingfile,".rds")
  if (!quiet){
    message("\n\nusing bigSNP file: ", bigsnp_path)
    message("with backing file: ", bigsnp_obj$genotypes$backingfile)
    message("make sure that you keep those files and don't delete them!")
  }

  indiv_meta <- as.list (indiv_meta)
  indiv_meta$genotypes <- new_vctrs_bigsnp(bigsnp_obj,
                                           bigsnp_file = bigsnp_path,
                                           indiv_id = bigsnp_obj$fam$sample.ID)

  tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )

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
                  bigsnp_file = bigsnp_file, # TODO is this redundant with the info in the bigSNP object?
                  bigsnp_file = tools::md5sum(bigsnp_file), # TODO make sure this does not take too long
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
    stop("loci does not include the compulsory columns 'name', 'chromosome', 'position','genetic_dist', allele_ref','allele_alt'")
  }
}


#' Create a bigSNP object from data.frames
#'
#' This function expects the indiv_meta and loci to have the correct columns
#' @param genotypes a genotype matrix
#' @param indiv_meta the individual meta information
#' @param loci the loci table
#' @keywords internal
gt_write_bigsnp_from_dfs <- function(genotypes, indiv_meta, loci,
                                      backingfile=NULL){

  if (is.null(backingfile)){
    backingfile <- tempfile()
  }
  check_valid_loci(loci)
  bigGeno <- bigstatsr::FBM.code256(
    nrow = nrow(genotypes),
    ncol = ncol(genotypes),
    code = bigsnpr::CODE_012,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )
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
  snp_list <- structure(list(genotypes = bigGeno, fam = fam, map = map),
                        class = "bigSNP")

  # save it and return the path of the saved object
  rds <- bigstatsr::sub_bk(bigGeno$backingfile, ".rds")
  saveRDS(snp_list, rds)
  return(snp_list)
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


##########################################
# convenient functs
.gt_bigsnp_cols <- function(.x){
  show_loci(.x)$big_index
}

.gt_bigsnp_rows <- function(.x){
  vctrs::vec_data(.x$genotypes)
}

