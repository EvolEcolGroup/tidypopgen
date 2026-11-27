# read in a vcf (with phased = TRUE, phased diploid genotypes are stored, see
# gen_tibble(); only with the cpp parser)
gen_tibble_vcf <- function(
  x,
  ...,
  parser = c("cpp", "vcfR"),
  n_cores = 1,
  chunk_size = NULL,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  backingfile = NULL,
  allow_duplicates = FALSE,
  names_as_int = FALSE,
  chromosomes_as_int = FALSE,
  quiet = FALSE,
  phased = FALSE
) {
  parser <- match.arg(parser)
  if (!file.exists(x)) {
    stop("x should be a valid file path pointing to a vcf: ", x)
  }

  if (parser == "cpp") {
    # check that the ellipses are empty (these are extra params for vcfR)
    if (length(list(...)) > 0) {
      stop("extra parameters can only be used with parser = 'vcfR'")
    }
    new_gen_tbl <- vcf_to_fbm_cpp(
      x,
      backingfile = backingfile,
      valid_alleles = valid_alleles,
      missing_alleles = missing_alleles,
      allow_duplicates = allow_duplicates,
      names_as_int = names_as_int,
      chromosomes_as_int = chromosomes_as_int,
      n_cores = n_cores,
      quiet = quiet,
      phased = phased
    )
  } else {
    new_gen_tbl <- vcf_to_fbm_vcfR(
      x,
      backingfile = backingfile,
      valid_alleles = valid_alleles,
      missing_alleles = missing_alleles,
      chunk_size = chunk_size,
      allow_duplicates = allow_duplicates,
      quiet = quiet,
      ...
    )
  }

  if (!quiet) {
    message("converting to a gen_tibble...")
  }
  return(new_gen_tbl)
}
