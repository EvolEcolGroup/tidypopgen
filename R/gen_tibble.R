#' Constructor for a `gen_tibble`
#'
#' A `gen_tibble` stores genotypes for individuals in a tidy format. DESCRIBE
#' here the format
#'
#' When loading packedancestry files, missing alleles will be converted from
#' 'X' to NA
#'
#' @param x can be:
#' - a string giving the path to a PLINK BED or PED file. The associated
#' BIM and FAM files for the BED, or MAP for PED are expected to be in the same
#' directory and have the same file name.
#' - a string giving the path to a RDS file storing a `bigSNP` object from
#' the `bigsnpr` package (usually created with [bigsnpr::snp_readBed()])
#' - a string giving the path to a vcf file. Note that we currently read the whole
#' vcf in memory with `vcfR`, so only smallish vcf can be imported. Only biallelic
#' SNPs will be considered.
#' - a string giving the path to a packedancestry .geno file. The associated
#' .ind and .snp files are expected to be in the same directory and share the
#' same file name prefix.
#' - a genotype matrix of dosages (0, 1, 2, NA) giving the dosage of the alternate
#' allele.
#' @param indiv_meta a list, data.frame or tibble with compulsory columns 'id'
#'  and 'population', plus any additional metadata of interest. This is only used
#' if `x` is a genotype matrix. Otherwise this information is extracted directly from
#' the files.
#' @param loci a data.frame or tibble, with compulsory columns 'name', 'chromosome',
#' and 'position','genetic_dist', 'allele_ref' and 'allele_alt'. This is only used
#' if `x` is a genotype matrix. Otherwise this information is extracted directly from
#' the files.
#' @param chunk_size the number of loci or individuals (depending on the format)
#'  processed at a time (currently used
#' if `x` is a vcf or packedancestry file)
#' @param ... if `x` is the name of a vcf file, additional arguments
#' passed to [vcfR::read.vcfR()]. Otherwise, unused.
#' @param parser the name of the parser used for VCF, either "cpp" to use
#' a fast C++ parser, or "vcfR" to use the R package `vcfR`. The latter is slower
#' but more robust; if "cpp" gives error, try using "vcfR" in case your VCF has
#' an unusual structure.
#' @param valid_alleles a vector of valid allele values; it defaults to 'A','T',
#' 'C' and 'G'.
#' @param missing_alleles a vector of values in the BIM file/loci dataframe that
#' indicate a missing value for the allele value (e.g. when we have a monomorphic
#' locus with only one allele). It defaults to '0' and '.' (the same as PLINK 1.9).
#' @param backingfile the path, including the file name without extension,
#' for backing files used to store the data (they will be given a .bk
#' and .RDS automatically). This is not needed if `x` is already an .RDS file.
#' If `x` is a .BED file and `backingfile` is left NULL, the backing file will
#' be saved in the same directory as the
#' bed file, using the same file name but with a different file type (.bk rather
#' than .bed). The same logic applies to .vcf files. If `x` is a genotype matrix and `backingfile` is NULL, then a
#' temporary file will be created (but note that R will delete it at the end of
#' the session!)
#' @param quiet provide information on the files used to store the data
#' @returns an object of the class `gen_tbl`.
#' @rdname gen_tibble
#' @export

gen_tibble <-
  function(x,
           ...,
           valid_alleles = c("A", "T", "C", "G"),
           missing_alleles = c("0","."),
           backingfile = NULL,
           quiet = FALSE) {

  UseMethod("gen_tibble", x)
}

###############################################################################
# character method for files
###############################################################################
#' @export
#' @rdname gen_tibble
gen_tibble.character <-
  function(x,
           ...,
           parser = c("vcfR","cpp"),
           chunk_size = NULL,
           valid_alleles = c("A", "T", "C", "G"),
           missing_alleles = c("0","."),
           backingfile = NULL,
           quiet = FALSE) {

    # parser for vcf
    parser <- match.arg(parser)


  # check that valid alleles does not contain zero
  if ("0" %in% valid_alleles){
    stop ("'0' can not be a valid allele (it is the default missing allele value!)")
  }

  if(is.null(backingfile)){
    backingfile <- change_duplicated_file_name(x)
  } else if(!is.null(backingfile)){
    backingfile <- change_duplicated_file_name(backingfile)
  }


  if ((tolower(file_ext(x))=="bed") || (tolower(file_ext(x))=="rds")){
    rlang::check_dots_empty()
    x_gt <- gen_tibble_bed_rds(x = x, ...,
                       valid_alleles= valid_alleles,
                       missing_alleles= missing_alleles,
                       backingfile = backingfile,
                       quiet = quiet)
  } else if ((tolower(file_ext(x))=="vcf") || (tolower(file_ext(x))=="gz")){
    return(gen_tibble_vcf(x = x, ..., parser = parser, chunk_size = chunk_size,
                   valid_alleles= valid_alleles,
                   missing_alleles= missing_alleles,
                   backingfile = backingfile, quiet = quiet))
  } else if (tolower(file_ext(x))=="ped"){
    x_gt <- gen_tibble_ped(x = x, ...,
                       valid_alleles= valid_alleles,
                       missing_alleles= missing_alleles,
                       backingfile = backingfile,
                       quiet = quiet)
  } else if (tolower(file_ext(x))=="geno"){
    x_gt <- gen_tibble_packedancestry(x = x, ...,
                   valid_alleles= valid_alleles,
                   missing_alleles= missing_alleles,
                   backingfile = backingfile,
                   quiet = quiet)
  } else  {
    stop("file_path should be pointing to a either a PLINK .bed or .ped file, a bigSNP .rds file or a VCF .vcf or .vcf.gz file")
  }

  loci <- show_loci(x_gt)
  new_loci <- add_chromosome_as_int(loci)
  show_loci(x_gt) <- new_loci

  file_in_use <- gt_save_light(x_gt, quiet = quiet)
  return(x_gt)
}


gen_tibble_bed_rds <- function(x, ...,
                               valid_alleles = c("A", "T", "C", "G"),
                               missing_alleles = c("0","."),
                               backingfile = NULL, quiet = FALSE){


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

  if (!("family.ID" %in% names(bigsnp_obj$fam))){
    indiv_meta <- list(id = bigsnp_obj$fam$sample.ID)
  } else {
    indiv_meta <- list(id = bigsnp_obj$fam$sample.ID,
                       population = bigsnp_obj$fam$family.ID)
  }
  # check if the bignsp_obj$fam table has ploidy column, if not, set ploidy to 2
  if ("ploidy" %in% names(bigsnp_obj$fam)){
    ploidy <- bigsnp_obj$fam$ploidy
  } else {
    ploidy <- 2
    bigsnp_obj$fam$ploidy <- 2
  }

  indiv_meta$genotypes <- new_vctrs_bigsnp(bigsnp_obj,
                                           bigsnp_file = bigsnp_path,
                                           indiv_id = bigsnp_obj$fam$sample.ID,
                                           ploidy = ploidy)

  # transfer some of the fam info to the metadata table if it is not missing (0 is the default missing value)
  fam_info <- .gt_get_bigsnp(indiv_meta)$fam
  if(!all(fam_info$paternal.ID==0)){
    indiv_meta$paternal_ID <- fam_info$paternal.ID
    indiv_meta$paternal_ID[indiv_meta$paternal_ID==0]<-NA
  }
  if(!all(fam_info$maternal.ID==0)){
    indiv_meta$maternal_ID <- fam_info$maternal.ID
    indiv_meta$maternal_ID[indiv_meta$maternal_ID==0]<-NA
  }
  # if sex is numeric
  if (inherits(fam_info$sex,"numeric")){
    if(!all(fam_info$sex==0)){
      indiv_meta$sex <-   dplyr::case_match(
        fam_info$sex,
        1 ~ "male",
        2 ~ "female",
        .default = NA,
        .ptype = factor(levels = c("female", "male"))
      )
    }
  }

  if (inherits(fam_info$affection,"numeric")){
  if(!all(fam_info$affection %in% c(0,-9))){
  indiv_meta$phenotype <- dplyr::case_match(
    fam_info$affection,
    1 ~ "control",
    2 ~ "case",
    -9 ~ NA,
    .default = NA,
    .ptype = factor(levels = c("control", "case"))
  )
  }
  }

  new_gen_tbl <- tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )
  check_allele_alphabet(new_gen_tbl, valid_alleles = valid_alleles,
                         missing_alleles = missing_alleles,
                         remove_on_fail = TRUE)
  show_loci(new_gen_tbl) <- harmonise_missing_values(show_loci(new_gen_tbl), missing_alleles = missing_alleles)
  return(new_gen_tbl)

}


###############################################################################
# matrix method to provide data directly from R
###############################################################################
#' @param ploidy the ploidy of the samples (either a single value, or
#' a vector of values for mixed ploidy). Only used if creating
#' a gen_tibble from a matrix of data; otherwise, ploidy is determined automatically
#' from the data as they are read.
#' @export
#' @rdname gen_tibble
gen_tibble.matrix <- function(x, indiv_meta, loci, ...,
                              ploidy = 2,
                              valid_alleles = c("A", "T", "C", "G"),
                              missing_alleles = c("0","."),
                              backingfile = NULL, quiet = FALSE){
  rlang::check_dots_empty()

  # check that valid alleles does not contain zero
  if ("0" %in% valid_alleles){
    stop ("'0' can not be a valid allele (it is the default missing allele value!)")
  }

  if (!inherits(loci, "data.frame") || inherits(x, "tbl")){
    stop("loci must be one of data.frame or tbl")
  }
  if (!inherits(indiv_meta, "data.frame") || inherits(x, "tbl") || is.list(x)){
    stop("indiv_meta must be one of data.frame, tbl, or list")
  }
  #if (!all(c("id", "population") %in% names(indiv_meta))){
  #  stop("ind_meta does not include the compulsory columns 'id' and 'population")
  #}
  if (!all(c("id") %in% names(indiv_meta))){
    stop("ind_meta does not include the compulsory column 'id")
  }
  # check that x (the genotypes) is numeric matrix
  if (inherits(x,"data.frame")){
    x <- as.matrix(x)
  }
  if (any(!inherits(x,"matrix"), !is.numeric(x))){
    stop("'x' is not a matrix of integers")
  }

  # check dimensions
  if (ncol(x)!=nrow(loci)){
    stop ("there is a mismatch between the number of loci in the genotype table x and in the loci table")
  }
  if (nrow(x)!=nrow(indiv_meta)){
    stop ("there is a mismatch between the number of loci in the genotype table x and in the loci table")
  }


  if(!is.null(backingfile)){
    backingfile <- change_duplicated_file_name(backingfile)
  }

  # use code for NA in FBM.256
#  x[is.na(x)]<-3

  bigsnp_obj <- gt_write_bigsnp_from_dfs(genotypes = x,
                                          indiv_meta = indiv_meta,
                                          loci = loci,
                                          backingfile = backingfile,
                                         ploidy = ploidy)

  bigsnp_path <- bigstatsr::sub_bk(bigsnp_obj$genotypes$backingfile,".rds")

  indiv_meta <- as.list (indiv_meta)
  indiv_meta$genotypes <- new_vctrs_bigsnp(bigsnp_obj,
                                           bigsnp_file = bigsnp_path,
                                           indiv_id = bigsnp_obj$fam$sample.ID,
                                           ploidy= ploidy)

  new_gen_tbl <- tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )
  check_allele_alphabet(new_gen_tbl, valid_alleles = valid_alleles,
                         missing_alleles = missing_alleles,
                        remove_on_fail = TRUE)
  show_loci(new_gen_tbl) <- harmonise_missing_values(show_loci(new_gen_tbl), missing_alleles = missing_alleles)

  loci <- show_loci(new_gen_tbl)
  new_loci <- add_chromosome_as_int(loci)
  show_loci(new_gen_tbl) <- new_loci

  files_in_use <- gt_save(new_gen_tbl, quiet = quiet)
  return(new_gen_tbl)

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
                                      backingfile=NULL,
                                      ploidy = ploidy){

  if (is.null(backingfile)){
    backingfile <- tempfile()
  }
  check_valid_loci(loci)
  # set up code (accounting for ploidy)
  code256 <- rep(NA_real_, 256)
  if (length(ploidy>1)){
    max_ploidy <- max(ploidy)
  } else {
    max_ploidy <- ploidy
  }

  if (is.na(max_ploidy)){
    stop("'ploidy' can not contain NAs")
  }

  code256[1:(max_ploidy+1)]<-seq(0,max_ploidy)

  # ensure max_ploidy is appropriate for the data
  if (any (genotypes>max_ploidy, na.rm=TRUE)){
    stop("max ploidy is set to ",max_ploidy," but genotypes contains indviduals with greater ploidy")
  }


  bigGeno <- bigstatsr::FBM.code256(
    nrow = nrow(genotypes),
    ncol = ncol(genotypes),
    code = code256,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )
  genotypes[is.na(genotypes)]<-max_ploidy+1
  bigGeno[] <- as.raw(genotypes)
  fam <- tibble(family.ID = indiv_meta$population,
                sample.ID = indiv_meta$id,
                paternal.ID = 0,
                maternal.ID = 0,
                sex = 0,
                affection = 0,
                ploidy = ploidy)

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

################################################################################
## vctrs_bigSNP class to store the genotype info in a gen_tibble
################################################################################

#' create a vctrs_bigSNP
#' @param bigsnp_obj the bigsnp object
#' @param bigsnp_file the file to which the bigsnp object was saved
#' @param indiv_id ids of individuals
#' @param ploidy the ploidy of the samples (either a single value, or
#' a vector of values for mixed ploidy).
#' @returns a vctrs_bigSNP object
#' @keywords internal
new_vctrs_bigsnp <- function(bigsnp_obj, bigsnp_file, indiv_id, ploidy = 2) {
  loci <- tibble::tibble(big_index = seq_len(nrow(bigsnp_obj$map)),
                         name = bigsnp_obj$map$marker.ID,
                         chromosome = bigsnp_obj$map$chromosome,
                         position = bigsnp_obj$map$physical.pos,
                         genetic_dist = bigsnp_obj$map$genetic.dist,
                         allele_ref = bigsnp_obj$map$allele2,
                         allele_alt = bigsnp_obj$map$allele1
  )

  if (length(unique(ploidy))>1){
    max_ploidy <- 0
  } else {
    max_ploidy <- max(ploidy)
  }
  vctrs::new_vctr(seq_len(nrow(bigsnp_obj$fam)),
                  bigsnp = bigsnp_obj,
                  bigsnp_file = bigsnp_file, # TODO is this redundant with the info in the bigSNP object?
                  bigsnp_md5sum = tools::md5sum(bigsnp_file), # TODO make sure this does not take too long
                  loci=loci,
                  names=indiv_id,
                  ploidy = max_ploidy,
                  class = "vctrs_bigSNP")
}

#' @export
summary.vctrs_bigSNP <- function(object, ...){
  summary(rep("bigSNP-genotypes",length(object)))
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
  if (!"genotypes" %in% names(.x)){
    stop("not a gen_tibble, 'genotype' column is missing")
  }
  if (!inherits(.x$genotypes,"vctrs_bigSNP")){
    stop("not a gen_tibble, the genotype column is not of class vctrs_bigSNP")
  }
}

# print method
#' @export
tbl_sum.gen_tbl <- function(x, ...) {
  c(
    "A gen_tibble" = paste(count_loci(x)," loci"),
    NextMethod()
  )
}


# function to check the allele alphabet
check_allele_alphabet <- function(x,
                                  valid_alleles = c("A", "T", "C", "G"),
                                  missing_alleles = c("0","."),
                                  remove_on_fail = FALSE){
  if (any(!show_loci(x)$allele_ref %in% c(valid_alleles,missing_alleles,NA),
          !show_loci(x)$allele_alt %in% c(valid_alleles,missing_alleles,NA))){
    if (remove_on_fail){ # remove files if they were generated
      if(file.exists(gt_get_file_names(x)[1])){
        file.remove(gt_get_file_names(x)[1])
      }
      if(file.exists(gt_get_file_names(x)[2])){
        file.remove(gt_get_file_names(x)[2])
      }
    }
    stop("valid alleles are ", paste(c(valid_alleles,missing_alleles), collapse=" ")," but ",
         paste(unique(c(show_loci(x)$allele_ref,show_loci(x)$allele_alt)), collapse=" "),
         " were found.")
  }

}

# set all missing values to NA
# loci_info is usually from show_loci()
harmonise_missing_values <- function (loci_info, missing_alleles =c("0",".")){
  # 0 is always considered as a missing value
  if (!"0" %in% missing_alleles){
    missing_alleles <- c(missing_alleles, "0")
  }
  loci_info$allele_ref[loci_info$allele_ref %in% missing_alleles]<-NA
  loci_info$allele_alt[loci_info$allele_alt %in% missing_alleles]<-NA
  return(loci_info)
}


# check for existing .bk files
change_duplicated_file_name <- function(file){

  file <- tools::file_path_sans_ext(file)

  bk <- paste0(file, ".bk")
  rds <- paste0(file, ".rds")

  if(file.exists(bk) && !file.exists(rds)){

    version <- 1

    base_name <- basename(file)

    version_pattern <- paste0(base_name, "_v(\\d+)\\.bk$")

    # read existing files to check for existing versions
  existing_files <- list.files(dirname(bk), pattern = paste0("^", base_name, "_v\\d+\\.bk$"))


    if (length(existing_files) > 0) {
      versions <- sub(version_pattern, "\\1", existing_files)
      versions <- as.numeric(versions)
      if (!any(is.na(versions))) {
        version <- max(versions) + 1
      }
    }

    new_file <- paste0(file,"_v",version)

    return(new_file)
  } else if (file.exists(bk) && file.exists(rds)){

    version <- 1

    base_name <- basename(file)

    version_pattern <- paste0(base_name, "_v(\\d+)\\.bk$")

    # read existing files to check for existing versions
    existing_files <- list.files(dirname(bk), pattern = paste0("^", base_name, "_v\\d+\\.bk$"))


    if (length(existing_files) > 0) {
      versions <- sub(version_pattern, "\\1", existing_files)
      versions <- as.numeric(versions)
      if (!any(is.na(versions))) {
        version <- max(versions) + 1
      }
    }

    new_file <- paste0(file,"_v",version)

    return(new_file)
  }

  return(file)

}

# adding a chr_int column
add_chromosome_as_int <- function(loci){


  if(is.integer(loci$chromosome)){
    loci$chr_int <- loci$chromosome

  } else if(is.numeric(loci$chromosome)){
    non_na_id <- !is.na(loci$chromosome)
    chr_int <- as.integer(loci$chromosome[non_na_id])

    loci$chr_int <- rep(NA_integer_, length(loci$chromosome))
    loci$chr_int[non_na_id] <- chr_int

  } else if(is.character(loci$chromosome)){
    non_na_id <- !is.na(loci$chromosome)
    chr_factor <- as.factor(loci$chromosome[non_na_id])
    chr_int <- as.integer(chr_factor)

    loci$chr_int <- rep(NA_integer_, length(loci$chromosome))
    loci$chr_int[non_na_id] <- chr_int


  } else if(is.factor(loci$chromosome)){
    non_na_id <- !is.na(loci$chromosome)
    chr_int <- as.integer(loci$chromosome[non_na_id])

    loci$chr_int <- rep(NA_integer_, length(loci$chromosome))
    loci$chr_int[non_na_id] <- chr_int

  } else {

    stop("Chromosome column should be integer, character, or factor")
  }
  return(loci)
}

