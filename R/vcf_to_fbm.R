#' Convert vcf to FBM.
#'
#' Convert a vcf file to a Filebacked Big Matrix (FBM) object.
#' This should work even for large vcf files that would not fit in memory.
#' TODO: this function is not yet complete.
#'
#' @param vcf_path the path to the vcf
#' @param chunks the chunk size to use on the vcf when loading the file
#' @param backingfile the name of the file to use as the backing file
#' @return path to the resulting rds file as class bigSNP.
#' @keywords internal

vcf_to_fbm <- function(
    vcf_path,
    chunks = 1000,
    backingfile = tempfile("vcf_matrix.bin")) {

  # count the variants in the file
  no_variants <- count_vcf_variants(vcf_path)
  # count individuals in the file
  no_individuals <- count_vcf_individuals(vcf_path)
  # set up chunks
  chunks_vec <- c(
    rep(chunks, floor(no_variants / chunks)),
    no_variants %% chunks
  )
  chunks_vec_index <- c(1, chunks_vec)

  # figure out ploidy from the first 1000 markers
  temp_vcf <- vcfR::read.vcfR(
    vcf_path,
    nrow = 1000
  )
  temp_gt <- vcfR::extract.gt(temp_vcf)
  ploidy <- apply(temp_gt,2,get_ploidy)
  max_ploidy <- max(ploidy, na.rm = TRUE) # remove NA in case we failed to figure out the ploidy of an individual

  # set up codes for the appropriate ploidy level
  code256 <- rep(NA_real_, 256)
  code256[1:(max_ploidy+1)]<-seq(0,max_ploidy)


  # create the file backed matrix
  file_backed_matrix <- bigstatsr::FBM.code256(
    nrow = no_individuals,
    ncol = 0,
    code = code256,
    backingfile = backingfile
  )

  for (i in seq_along(chunks_vec)) {
    temp_vcf <- vcfR::read.vcfR(
      vcf_path,
      nrow = chunks_vec[i],
      skip = sum(chunks_vec[seq_len(i-1)])
    )
    # filter any marker that is not biallelic
    bi <- vcfR::is.biallelic(temp_vcf)
    gt <- vcfR::extract.gt(temp_vcf)
    gt <- gt[bi,,drop=FALSE]
    if (nrow(gt)>1){
      gt <- t(apply(gt,2,poly_indiv_dosage, max_ploidy=max_ploidy))
    } else if (nrow(gt)==1){ # if we only have one marker
      gt <- matrix(apply(gt,2,poly_indiv_dosage, max_ploidy=max_ploidy),ncol=1)
    } else  {
      next
    }
    # expand the file backed matrix according to the size of the gt matrix
    # get current size
    index_start <- dim(file_backed_matrix)[2]+1
    # add the new columns
    file_backed_matrix$add_columns(ncol(gt))
    # fill them in
    file_backed_matrix[
      ,
      index_start:(index_start+ncol(gt)-1)
    ] <- gt
  }

  bigsnp_save <- structure(list(
    genotypes = file_backed_matrix,
    # TODO: these will have to be filled in?
    # YES they do! copy what we do in gen_tibble_vcf
    fam = NULL,
    map = NULL
  ), class = "bigSNP")

  # add .rds extension to backingfile
  out <- paste0(backingfile, ".rds")
  saveRDS(bigsnp_save, file = out)
  # and return the path to the rds
  out
}

# get ploidy for a given individual
get_ploidy <- function(x){
  max(sapply(strsplit(x,"[/|]"),function(x) length(x) ))
}

# get dosages for all the genotypes of an individual (x is a vector of genotypes in the standard vcf format)
poly_indiv_dosage <- function (x, max_ploidy){
  sapply(strsplit(x,"[/|]"),poly_genotype_dosage , max_ploidy)
}

# get dosages for a genotype x (as a vector of 0 and 1)
poly_genotype_dosage <- function (x, max_ploidy){
  if (!is.na(x[1]) && x[1]!="."){
    as.raw(sum(as.numeric(x))+1)
  } else {
    return(as.raw(max_ploidy+2))
  }
}

