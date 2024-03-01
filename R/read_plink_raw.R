#' Read data in PLINK raw format and store it to a `gen_tibble`
#'
#' This function reads a PLINK raw file (and optionally a map file) and stores
#' the data as a `gen_tibble`. The raw file should be generated with the flags
#' `--include-alt` and `--recodeA`. The map file is generated with `--recode` (in other words,
#' you will need to run PLINK twice to generate the necessary files).
#' @param file a character string giving the path to the file to convert,
#' with the extension ".raw"
#' @param map_file character string indicating the path to a ".map" file,
#' which contains information about the SNPs (chromosome, position).
#' @param quiet logical stating whether conversion messages should be printed
#'  (TRUE, default) or not (FALSE).
#' @param chunk_size an integer indicating the number of genomes to be read at a time;
#' larger values require more RAM but decrease the time needed to read the data.
#' @param parallel a logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE, default), or not (FALSE);
#' requires the package parallel to be installed (see details).
#' @param n_cores if parallel is TRUE, the number of cores to be used in the
#' computations; if NULL, then the maximum number of cores
#' available on the computer is used.
#' @param ... other arguments to be passed to other functions - currently not used.
#' @returns an object of the class `gen_tbl`.
#' @examples
#' raw_path <- system.file("extdata/san.raw", package = "tidypopgen")
#' map_path <- system.file("extdata/san.map", package = "tidypopgen")
#' san_gen_tbl <- read_plink_raw(file = raw_path, map_file = map_path)
#'
#' @export


read_plink_raw <- function(file, map_file=NULL, quiet=FALSE, chunk_size=1000,
                       parallel=FALSE, n_cores=NULL, ...){

  ext <- adegenet::.readExt(file)
  ext <- toupper(ext)
  if(ext != "RAW") warning("wrong file extension - '.raw' expected")
  if(!quiet) cat("Reading PLINK raw format into a gen_tibble... \n\n")
  if(parallel && is.null(n.cores)){
    n.cores <- parallel::detectCores()
  }

  if(!quiet) cat("Reading loci information... \n")
  col_names <- scan(file,what="character",sep=" ",quiet=TRUE,  nlines=1, blank.lines.skip=FALSE)
  ind_meta_list <- lapply(1:6,function(i) NULL)
  names(ind_meta_list) <- col_names[1:6]
  loci_names <- col_names[7:length(col_names)]
  # remove underscore followed by a digit at the end of locus name
  # when would this happen?!?
#  loci_names <- gsub("_[1-9]$","",loci_names)

  # now parse the alleles
  #TODO we should check that we have info for both alleles

  # parse the information from the loci col names
  # the format is rs3131972_A(/G)
  loci <- tibble::tibble(name = substr(loci_names,1,nchar(loci_names)-6),
                 chromosome = NA, position = NA,
                 allele_ref = tolower(substr(loci_names,nchar(loci_names)-1,nchar(loci_names)-1)),
                 allele_alt = tolower(substr(loci_names,nchar(loci_names)-4,nchar(loci_names)-4)))
  # replace zeroes (missing) with NA
  loci$allele_alt[loci$allele_alt=="0"] <- NA

  if (!is.null(map_file)){
    loci_map <- utils::read.table(map_file)
    names(loci_map) <- c("chromosome","name","cM","position")
    if(!all(loci$name==loci_map$name)){
      stop("there is a mismatch between the raw and the map file, the loci are not the same")
    }
    loci$chromosome <- factor(loci_map$chromosome)
    loci$position <- loci_map$position
    loci$cM <- loci_map$cM #TODO should we turn 0 into missing??
  }




  if(!quiet) cat("Reading and converting genotypes... \n")

  res <- list() # this will be a list of SNPbin objects

  ## initialize reading
  lines.to.skip <- 1
  txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=chunk_size)
  txt <- lapply(txt, function(e) unlist(strsplit(e,"[[:blank:]]+") ))

  COUNT <- 0 # used to count the nb reads

  while(length(txt)>0){
    COUNT <- COUNT + 1
    if(!quiet) {
      if(COUNT %% 5 == 0){
        cat(length(res)+length(txt))
      } else {
        cat(".")
      }
    }


    ## handle misc info
    temp <- lapply(txt, function(e) e[1:6])
    for(i in 1:6){
      ind_meta_list[[i]] <- c(ind_meta_list[[i]], unlist(lapply(temp, function(e) e[[i]])) )
    }


    ## build SNPbin objects
    txt <- lapply(txt, function(e) suppressWarnings(as.integer(e[-(1:6)])))

    if(parallel){
      res <- c(res, parallel::mclapply(txt, function(e) methods::new("SNPbin", snp=e, ploidy=2L),
                                       mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE) )
    } else {
      res <- c(res, lapply(txt, function(e) methods::new("SNPbin", snp=e, ploidy=2L)) )
    }

    lines.to.skip <-lines.to.skip + length(txt)

    ## read lines
    txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=chunk_size)
    txt <- lapply(txt, function(e) unlist(strsplit(e,"[[:blank:]]+") ))
  }


  ## MAKE A FEW CHECKS ##
  if(!all(sapply(res, adegenet::nLoc)==nrow(loci))) {
    stop(paste("some individuals do not have",nrow(loci),"SNPs."))
  }


  if(!quiet) cat("Reading individual metadata... \n")
  names(ind_meta_list) <- c("population","id","pat","mat","sex","phenotype")
  ind_meta_list <- ind_meta_list[c("id","population","sex","pat","mat","phenotype")]
  # recode some of these values
  ind_meta_list$population <- factor(ind_meta_list$population)
  ind_meta_list$pat <- as.character(dplyr::recode(ind_meta_list$pat,"0" = NA))
  ind_meta_list$mat <- as.character(dplyr::recode(ind_meta_list$mat,"0" = NA))
  ind_meta_list$phenotype <- dplyr::case_match(
    ind_meta_list$phenotype,
    "1" ~ "control",
    "2" ~ "case",
    "-9" ~ NA,
    .default = NA,
    .ptype = factor(levels = c("control", "case"))
  )
  ind_meta_list$sex <-   dplyr::case_match(
    ind_meta_list$sex,
    "1" ~ "male",
    "2" ~ "female",
    .default = NA,
    .ptype = factor(levels = c("female", "male"))
  )
  ind_meta_list$genotypes <- res
  attr(ind_meta_list$genotypes,"loci")<-tibble::as_tibble(loci)

  if(!quiet) cat("Building final object... \n")
  res <- tibble::new_tibble(
    ind_meta_list,
    class = "gen_tbl"
  )


  if(!quiet) cat("...done.\n")
  return(res)
}

