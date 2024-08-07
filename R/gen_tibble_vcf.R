# read in a vcf
gen_tibble_vcf <- function(x, ..., parser = "vcfR",
                           chunk_size = NULL,
                           valid_alleles = c("A", "T", "C", "G"),
                           missing_alleles = c("0","."),
                           backingfile = NULL, quiet = FALSE) {
  if (parser == "cpp"){
    rds_path <- vcf_to_fbm_cpp(x,
                                backingfile = backingfile,
                                chunk_size = chunk_size,
                                quiet = quiet)
  } else {
    rds_path <- vcf_to_fbm_vcfR(x,
                                backingfile = backingfile,
                                chunk_size = chunk_size,
                                quiet = quiet)
  }

  if (!quiet){
    message("converting to a gen_tibble...")
  }
  gen_tibble(rds_path,
             valid_alleles = valid_alleles,
             missing_alleles = missing_alleles,
             backingfile = backingfile,
             quiet = quiet)
}
