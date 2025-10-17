# A function to read geno packedancestrymap files
gen_tibble_packedancestry <- function(
    x,
    ...,
    valid_alleles = c("A", "T", "C", "G"),
    missing_alleles = c("0", "."),
    allow_duplicates = FALSE,
    backingfile = NULL,
    quiet = FALSE) {
  # Substitute .geno with .snp
  map_file <- sub("\\.geno$", ".snp", x)
  if (!file.exists(map_file)) {
    stop("snp file ", map_file, " does not exist")
  }
  ind_file <- sub("\\.geno$", ".ind", x)
  if (!file.exists(ind_file)) {
    stop("ind file ", ind_file, " does not exist")
  }

  if (is.null(backingfile)) {
    backingfile <- sub("\\.geno$", "", x)
  }
  if (file_ext(backingfile) == "bk") {
    backingfile <- bigstatsr::sub_bk(backingfile, "")
  }

  # read the packed ancestry in chunks
  # count the individuals
  indiv_table <- utils::read.table(ind_file, header = FALSE)
  names(indiv_table) <- c("id", "sex", "population")
  no_individuals <- nrow(indiv_table)

  indiv_meta <- list(
    id = indiv_table$id,
    population = indiv_table$population
  )

  indiv_meta <- validate_indiv_meta(as.data.frame(indiv_meta))

  # count the loci
  loci_table <- utils::read.table(map_file, header = FALSE)
  names(loci_table) <- c(
    "name",
    "chromosome",
    "genetic_dist",
    "position",
    "allele_ref",
    "allele_alt"
  )
  loci_table$allele_alt[loci_table$allele_alt == "X"] <- NA
  no_variants <- nrow(loci_table)

  loci <- validate_loci(loci_table,
    check_alphabet = TRUE,
    harmonise_loci = TRUE,
    check_duplicates = TRUE,
    allow_duplicates = allow_duplicates,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles
  )

  # verify that these numbers are compatible with the geno file
  conn <- file(x, "rb")
  hd <- strsplit(readBin(conn, "character", n = 1), " +")[[1]]
  close(conn)
  if (no_individuals != as.numeric(hd[2])) {
    stop(
      "Number of individuals in geno file does not match the number of ",
      "individuals in the ind file"
    )
  }
  if (no_variants != as.numeric(hd[3])) {
    stop(
      "Number of variants in geno file does not match the number of ",
      "variants in the snp file"
    )
  }

  # create a matrix to store the data
  file_backed_matrix <- bigstatsr::FBM.code256(
    nrow = no_individuals,
    ncol = no_variants,
    code = bigsnpr::CODE_012,
    backingfile = backingfile
  )

  # Fill the FBM from bedfile
  reach_eof <- read_packedancestry(x, file_backed_matrix,
    tab = get_packedancestry_code()
  )

  if (!reach_eof) warning("EOF of bedfile has not been reached.")

  # construct path
  fbm_path <- bigstatsr::sub_bk(file_backed_matrix$backingfile, ".rds")

  # create genotypes column
  indiv_meta$genotypes <- new_vctrs_bigsnp(
    fbm_obj = file_backed_matrix,
    fbm_file = fbm_path,
    loci = loci,
    indiv_id = indiv_meta$id,
    ploidy = 2,
    fbm_ploidy = rep(2, nrow(indiv_meta))
  )

  new_gen_tbl <- tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )
  return(new_gen_tbl)
}
