#' Export a `gen_tibble` object to PLINK ped format
#'
#' This function exports all the information of a `gen_tibble` object into
#' a PLINK ped file.
#'
#' This function generates a PLINK raw file with alternate allele information
#' (as outputted using `--include-alt` and
#' `--recodeA`), as well as a map file (as generated with `--recode`).
#'
#' @param x a [`gen_tibble`] object
#' @param file a character string giving the path to the file to convert,
#' with the extension ".raw" or "ped"
#' @param plink_format character, one of "raw" or "ped"
#' @param chunk_size very large objects are written in chunks to minimise the
#' memory footprint. Increase `chunk_site` for higher speed, decrease it for
#' smaller memory footprint (defaults to 10000)
#' @param overwrite boolean whether to overwrite the file
#' @returns TRUE if successful
#' @export


gt_as_plink <- function(x, file, plink_format = c("raw","ped"), chunk_size = 10000, overwrite = TRUE){

  plink_format <- match.arg(plink_format)
  if (tolower(adegenet::.readExt(file))!=plink_format){
    file <- paste0(file,".",plink_format)
  }

  if (file.exists(file)){
    if (overwrite){
      file.remove(file)
    } else {
      stop("file already exists; remove if first or set 'overwrite' = TRUE")
    }
  }

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
  return(TRUE)
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
  dplyr::case_match(
    x,
    "0" ~ genotypes[1],
    "1" ~ genotypes[2],
    "2" ~ genotypes[3]
  )

}
