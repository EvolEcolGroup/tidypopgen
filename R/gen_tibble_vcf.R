# read in a vcf
gen_tibble_vcf <- function(x, ...,
                           loci_per_chunk = 100000,
                           valid_alleles = c("A", "T", "C", "G"),
                           missing_alleles = c("0","."),
                           backingfile = NULL, quiet = FALSE) {
  rds_path <- vcf_to_fbm(x,backingfile = backingfile,loci_per_chunk = loci_per_chunk, quiet = quiet)
  gen_tibble(rds_path, valid_alleles = valid_alleles, missing_alleles = missing_alleles,
             backingfile = backingfile, quiet = quiet)
}
