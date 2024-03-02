#' Export a `gen_tibble` object to PLINK raw format
#'
#' This function exports all the information of a `gen_tibble` object into
#' a PLINK raw file.
#'
#' This function generates a PLINK raw file with alternate allele information
#' (as outputted using `--include-alt` and
#' `--recodeA`), as well as a map file (as generated with `--recode`).
#'
#' @param x a [`gen_tibble`] object
#' @param file a character string giving the path to the file to convert,
#' with the extension ".raw"
#' @param chunk_size very large objects are written in chunks to minimise the
#' memory footprint. Increase `chunk_site` for higher speed, decrease it for
#' smaller memory footprint (defaults to 10000)
#' @param overwrite boolean whether to overwrite the file
#' @returns TRUE if successful
#' @export


write_plink_raw <- function(x, file, chunk_size = 10000, overwrite = TRUE){

  if (toupper(adegenet::.readExt(file))!="RAW"){
    file <- paste0(file,".raw")
  }
  if (file.exists(file)){
    if (overwrite){
      file.remove(file)
    } else {
      stop("file already exists; remove if first or set 'overwrite' = TRUE")
    }

  }

  # format the col names
  c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  loci <- show_loci(x)
  # replace missing value with zero
  loci$allele_alt[is.na(loci$allele_alt)]<-"0"
  # create col names to use in the file
  raw_col_names<- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE",
  paste0(loci$name,"_",toupper(loci$allele_alt),"(/",toupper(loci$allele_ref),")"))

  n_ind <- nrow(x)

  chunks <- split(1:n_ind, ceiling(seq_along(1:n_ind)/chunk_size))

  # create the fam info
  ind_meta <- tibble(population = x$population,
                     id = x$id,
                     pat = pull_NA(x, "pat"),
                     mat = pull_NA(x, "mat"),
                     sex = pull_NA(x, "sex"),
                     phenotype = pull_NA(x, "phenotype")
                     )
  # recode some variables
  ind_meta$sex<- dplyr::case_match(as.character(ind_meta$sex),'female'~'2','male'~'1',.default="0")
  ind_meta$pat[is.na(ind_meta$pat)]<-0
  ind_meta$mat[is.na(ind_meta$mat)]<-0
  ind_meta$phenotype<-dplyr::case_match(as.character(ind_meta$phenotype),
                    'control'~"1",'case'~"2", .default=ind_meta$phenotype )
  ind_meta$phenotype[is.na(ind_meta$phenotype)]<--9
  # loop over to write
  for (i_chunk in seq_len(length(chunks))){
    chunk <- chunks[[i_chunk]]
    raw_table <- cbind(ind_meta[chunk,],
    show_genotypes(x[chunk,]))
    colnames(raw_table) <- raw_col_names
    # append column names only the first time, when the file does not exist
    utils::write.table(raw_table,
                file=file,
                sep = " ",
                row.names = FALSE,
                col.names=!file.exists(file),
                append = file.exists(file),
                quote = FALSE)

  }
  ## now write the map file

  loci_meta <- show_loci(x)
  map_file <- paste0(substr(file, 1, nchar(file)-4),".map")
  map_table <- data.frame(chrom = pull_NA(x,"chromosome"),
                          id = pull_NA(x,"name"),
                          cM = pull_NA(x,"cM"),
                          position = pull_NA(x,"position"))
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
