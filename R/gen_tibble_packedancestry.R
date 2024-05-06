#A function to read geno packedancestry files
gen_tibble_geno <- function(x, ...,
                            valid_alleles = c("A", "T", "C", "G"),
                            missing_alleles = c("0","."),
                            backingfile = NULL, quiet = FALSE) {
  # Substitute .ped with .map
  map_file <- sub("\\.geno$", ".snp", x)
  if (!file.exists(map_file)){
    stop("snp file ",map_file," does not exist")
  }
  ind_file <- sub("\\.geno$", ".ind", x)
  if (!file.exists(ind_file)){
    stop("ind file ",ind_file," does not exist")
  }

  res <- admixtools::read_packedancestrymap(sub("\\.geno$", "", x),
                                                  transpose = TRUE,
                                            ...)
  names(res$ind)<-c("id","sex","population")
  #TODO check that allele_ref and allele_alt are not swapped
  names(res$snp)<-c("name", "chromosome",'genetic_dist','position', 'allele_ref','allele_alt')


  # using the gen_tibble.matrix method
  new_gen_tbl <- gen_tibble(x = res$geno,
                            indiv_meta = res$ind,
                            loci = res$snp,
                            backingfile = backingfile, quiet=quiet)
  check_allele_alphabet (new_gen_tbl, valid_alleles = valid_alleles,
                         missing_alleles = missing_alleles)
  show_loci(new_gen_tbl) <- harmonise_missing_values(show_loci(new_gen_tbl), missing_alleles = missing_alleles)
  return(new_gen_tbl)

}
