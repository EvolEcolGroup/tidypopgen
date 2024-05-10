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
  no_individuals <- count_vcf_individuals(vcf_path)
  chunks_vec <- c(
    rep(chunks, floor(no_variants / chunks)),
    no_variants %% chunks
  )
  chunks_vec_index <- c(1, chunks_vec)

  # create the file backed matrix
  file_backed_matrix <- bigstatsr::FBM.code256(
    nrow = no_individuals,
    ncol = 0,
    # TODO: when ploidy code comes in, this will have to be changed
    code = bigsnpr::CODE_012,
    backingfile = backingfile
  )

  for (i in seq_along(chunks_vec)) {
    temp_vcf <- vcfR::read.vcfR(
      vcf_path,
      nrow = chunks_vec[i],
      skip = sum(chunks_vec[1:(i - 1)])
    )

    gt <- vcfR::extract.gt(temp_vcf)
    gt[gt == "0|0"] <- 0
    gt[gt == "0|1"] <- 1
    gt[gt == "1|0"] <- 1
    gt[gt == "1|1"] <- 2
    gt[gt == "0/0"] <- 0
    gt[gt == "0/1"] <- 1
    gt[gt == "1/0"] <- 1
    gt[gt == "1/1"] <- 2
    # FIXME: this is currently generating loads of warnings.
    gt <- apply(gt, 2, as.raw)

    file_backed_matrix$add_columns(chunks_vec[i])
    file_backed_matrix[
      ,
      cumsum(chunks_vec_index)[i]:cumsum(chunks_vec)[i]
    ] <- t(gt)
  }

  bigsnp_save <- structure(list(
    genotypes = file_backed_matrix,
    # TODO: these will have to be filled in?
    fam = NULL,
    map = NULL
  ), class = "bigSNP")

  # add .rds extension to backingfile
  out <- paste0(backingfile, ".rds")
  saveRDS(bigsnp_save, file = out)
  # and return the path to the rds
  out
}
