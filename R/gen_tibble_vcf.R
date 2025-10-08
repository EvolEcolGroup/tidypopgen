# read in a vcf
gen_tibble_vcf <- function(
    x,
    ...,
    parser = c("cpp", "vcfR"),
    n_cores = 1, # ignored by this function, there is no multithreading
    chunk_size = NULL,
    valid_alleles = c("A", "T", "C", "G"),
    missing_alleles = c("0", "."),
    backingfile = NULL,
    allow_duplicates = FALSE,
    quiet = FALSE) {
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
      quiet = quiet
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
