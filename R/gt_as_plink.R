#' Export a `gen_tibble` object to PLINK bed format
#'
#' This function exports all the information of a `gen_tibble` object into
#' a PLINK bed, ped or raw file (and associated files, i.e. .bim and .fam for .bed;
#' .fam for .ped).
#'
#' @param x a [`gen_tibble`] object
#' @param file a character string giving the path to output file. If left to NULL,
#' the output file will have the same path and prefix of the backingfile.
#' @param type one of "bed", "ped" or "raw"
#' @param overwrite boolean whether to overwrite the file.
#' @returns the path of the saved file
#' @export


gt_as_plink <- function(x, file = NULL, type = c("bed","ped","raw"),
                           overwrite = TRUE){
  type <- match.arg(type)

  if (is.null(file)){
    file <- bigstatsr::sub_bk(attr(x$genotypes,"bigsnp")$genotypes$backingfile,paste0(".",type))
  }
  if (file_ext(file)!=type){
    file <- paste0(file,".",type)
  }

 if (type=="bed"){
    all_files <- c(file,
                   bigsnpr::sub_bed(file,".bim"),
                   bigsnpr::sub_bed(file,".fam")
                   )
  } else if (type=="ped"){
    all_files <- c(file,
                   gsub(".ped",".map",file))
  } else if (type=="raw"){
    all_files <- file
  }

  if (any(file.exists(all_files))){
    if (overwrite){
      file.remove(all_files[file.exists(all_files)])
  } else {
      stop("at least one of", all_files," already exists; remove if first or set 'overwrite' = TRUE")
    }
  }

  if (type == "bed"){
    gt_write_bed(x, file)
  } else {
    gt_write_ped_raw(x,file,type)
  }
}


gt_write_bed <- function(x, file) {
  bed_path <- bigsnpr::snp_writeBed(attr(x$genotypes,"bigsnp"),
                        bedfile = file,
                        ind.row = vctrs::vec_data(x$genotypes),
                        ind.col = show_loci(x)$big_index)
  # the bim and fam file only contain the original information in the bigSNP object
  # TODO we should update them with the info from the gentibble
  bed_path
}



gt_write_ped_raw <- function(x, file, plink_format = c("raw","ped"), chunk_size = 10000){

  # loci information
  loci <- show_loci(x)
  # replace missing value with zero
  loci$allele_alt[is.na(loci$allele_alt)]<-"0"

  # create col names to use in the raw file
  raw_col_names<- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE",
                    paste0(loci$name,"_",toupper(loci$allele_alt),"(/",toupper(loci$allele_ref),")"))

  # create the info for the fam file
  indiv_meta <- tibble(population = x$population,
                       id = x$id,
                       pat = pull_NA(x, "pat"),
                       mat = pull_NA(x, "mat"),
                       sex = pull_NA(x, "sex"),
                       phenotype = pull_NA(x, "phenotype")
  )
  # recode some variables
  indiv_meta$sex<- dplyr::case_match(as.character(indiv_meta$sex),'female'~'2','male'~'1',.default="0")
  indiv_meta$pat[is.na(indiv_meta$pat)]<-0
  indiv_meta$mat[is.na(indiv_meta$mat)]<-0
  indiv_meta$phenotype<-dplyr::case_match(as.character(indiv_meta$phenotype),
                                          'control'~"1",'case'~"2", .default=indiv_meta$phenotype )
  indiv_meta$phenotype[is.na(indiv_meta$phenotype)]<--9

  # create chunks
  n_ind <- nrow(x)
  chunks <- split(1:n_ind, ceiling(seq_along(1:n_ind)/chunk_size))

  # loop over to write
  for (i_chunk in seq_len(length(chunks))){
    chunk <- chunks[[i_chunk]]
    raw_table <- cbind(indiv_meta[chunk,],
                       show_genotypes(x[chunk,]))

    # now recode the genotypes with letters if raw
    if (plink_format=="ped"){
      for (i in 1:(ncol(raw_table)-6)){
        raw_table[,i+6]<-recode_genotype(raw_table[,i+6], loci$allele_ref[i], loci$allele_alt[i])
      }
    }

    colnames(raw_table) <- raw_col_names
    # append column names only the first time, when the file does not exist
    utils::write.table(raw_table,
                       file=file,
                       sep = " ",
                       row.names = FALSE,
                       col.names=(!file.exists(file) & plink_format=="raw"),
                       append = file.exists(file),
                       quote = FALSE)

  }

  ## create info for the map file
  loci_meta <- show_loci(x)
  map_file <- paste0(substr(file, 1, nchar(file)-4),".map")
  map_table <- tibble(chrom = pull_NA(loci_meta,"chromosome"),
                      id = pull_NA(loci_meta,"name"),
                      cM = pull_NA(loci_meta,"cM"),
                      position = pull_NA(loci_meta,"position"))
  map_table[is.na(map_table)]<-0

  utils::write.table(
    map_table,
    file = map_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  return(file)
}


# internal function returning either a vector from a given column or NA if it does not exist
pull_NA <- function(.x, .col_name) {
  if (.col_name %in% names(.x)){
    return(.x %>% pull(.col_name))
  } else {
    return(rep(NA,nrow(.x)))
  }
}

# recode dosage as letter genotypes
recode_genotype <- function(x, allele_ref, allele_alt){
  x <-as.character(x)
  genotypes <- c(paste(allele_ref, allele_ref),
                 paste(allele_ref, allele_alt),
                 paste(allele_alt, allele_alt))
  x <- dplyr::case_match(
    x,
    "0" ~ genotypes[1],
    "1" ~ genotypes[2],
    "2" ~ genotypes[3]
  )
  x[is.na(x)]<-c("0 0")
  x

}


