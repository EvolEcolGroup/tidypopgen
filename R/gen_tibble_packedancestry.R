# A function to read packedancestrymap files (.geno, .snp and .ind), with
# genotypes in either of the two layouts of EIGENSOFT: GENO (SNP-major) or
# TGENO (individual-major, as in recent releases of the AADR). The layout is
# read from the header of the genotype file, whatever its extension (.geno or
# .tgeno)
gen_tibble_packedancestry <- function(
  x,
  ...,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  allow_duplicates = FALSE,
  backingfile = NULL,
  n_cores = 1,
  quiet = FALSE
) {
  geno_ext <- file_ext(x)
  # Substitute .geno (or .tgeno) with .snp and .ind
  map_file <- sub_ext_ignore_case(x, geno_ext, ".snp")
  if (!file.exists(map_file)) {
    stop("snp file ", map_file, " does not exist")
  }
  ind_file <- sub_ext_ignore_case(x, geno_ext, ".ind")
  if (!file.exists(ind_file)) {
    stop("ind file ", ind_file, " does not exist")
  }

  if (is.null(backingfile)) {
    backingfile <- sub("\\.t?geno$", "", x, ignore.case = TRUE)
  }
  if (file_ext(backingfile) == "bk") {
    backingfile <- bigstatsr::sub_bk(backingfile, "")
  }

  # check the header of the genotype file before reading the other files
  geno_header <- read_packedancestry_header(x)

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
  # (names and alleles are read as characters)
  loci_table <- read_loci_file(map_file, character_cols = c(1, 5, 6))
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

  # read the genotypes into a new FBM (the object that wrote it is kept, as a
  # new one from its .rds file would map the backing file again, and
  # unmapping the first one takes seconds for large files)
  file_backed_matrix <- fbm_read_packedancestry(x,
    n_indiv = no_individuals,
    n_snp = no_variants,
    backingfile = backingfile,
    n_cores = n_cores,
    header = geno_header
  )
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

# Read the genotypes of a packedancestrymap file (GENO or TGENO) into a new
# FBM.code256 (with the code of bigsnpr, 0, 1, 2 and NA, giving the dosage of
# the second allele of the .snp file), saved next to its backing file. The
# file is checked (with check_packedancestry_file()) before the backing file is
# created, so that a bad file leaves nothing behind. The genotypes are decoded
# (GENO) or transposed (TGENO) in parallel with n_cores cores
# (read_packedancestry_to_fbm(), in src/read_packedancestry.cpp). For a TGENO
# file, the whole file is read into memory if it is smaller than
# max_buffer_bytes, and otherwise it is read in slabs of loci, using up to
# max_buffer_bytes of memory. Returns the FBM.
fbm_read_packedancestry <- function(geno_file,
                                    n_indiv,
                                    n_snp,
                                    backingfile,
                                    n_cores = 1,
                                    header = NULL,
                                    max_buffer_bytes = 2^31) {
  backingfile <- path.expand(backingfile)
  geno_file <- path.expand(geno_file)
  if (file.exists(paste0(backingfile, ".bk"))) {
    stop("the backing file ", paste0(backingfile, ".bk"), " already exists")
  }
  if (file.exists(paste0(backingfile, ".rds"))) {
    stop("the backing file ", paste0(backingfile, ".rds"), " already exists")
  }
  geno_format <- check_packedancestry_file(geno_file,
    n_indiv = n_indiv,
    n_snp = n_snp,
    header = header
  )
  big_geno <- bigstatsr::FBM.code256(
    nrow = n_indiv,
    ncol = n_snp,
    code = bigsnpr::CODE_012,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )
  read_packedancestry_to_fbm(
    geno_file,
    big_geno,
    get_packedancestry_code(),
    transposed = (geno_format == "TGENO"),
    n_threads = n_cores,
    max_buffer_bytes = max_buffer_bytes
  )
  saveRDS(big_geno, bigstatsr::sub_bk(big_geno$backingfile, ".rds"))
  big_geno
}

# Read the header of a packedancestrymap genotype file: a string such as
# "GENO 5 16 be181b4d 9f75a306" or "TGENO 5 16 be181b4d 9f75a306" (the
# layout, the numbers of individuals and of loci, and hashes of their names),
# padded with NUL bytes to at least 48 bytes. Returns a list with the format
# ("GENO", SNP-major, or "TGENO", individual-major), n_indiv and n_snp.
read_packedancestry_header <- function(geno_file) {
  bytes <- readBin(geno_file, what = "raw", n = 48)
  is_magic <- function(magic) {
    magic <- charToRaw(magic)
    length(bytes) >= length(magic) &&
      identical(bytes[seq_along(magic)], magic)
  }
  # (TGENO first, as its magic does not start with GENO)
  if (is_magic("TGENO ")) {
    geno_format <- "TGENO"
  } else if (is_magic("GENO ")) {
    geno_format <- "GENO"
  } else {
    # EIGENSTRAT (and LEA) .geno files are text, with a digit (0, 1, 2 or 9)
    # for each genotype
    if (length(bytes) > 0 &&
          all(bytes %in% charToRaw("0129\r\n"))) {
      stop(
        "the file ", geno_file, " looks like an EIGENSTRAT (text) .geno ",
        "file, which is not supported; convert it to packedancestrymap ",
        "(e.g. with convertf from EIGENSOFT, using 'outputformat: ",
        "PACKEDANCESTRYMAP')"
      )
    }
    stop(
      "the file ", geno_file, " is not a packedancestrymap genotype file ",
      "(its header should start with GENO or TGENO)"
    )
  }
  # the header string ends at the first NUL byte
  nul <- which(bytes == as.raw(0))
  if (length(nul) > 0) {
    bytes <- bytes[seq_len(nul[1] - 1)]
  }
  fields <- strsplit(rawToChar(bytes), "[[:space:]]+")[[1]]
  counts <- suppressWarnings(as.numeric(fields[2:3]))
  if (length(fields) < 3 || anyNA(counts) || any(counts < 0) ||
        any(counts != round(counts))) {
    stop(
      "the header of the ", geno_format, " file ", geno_file,
      " is malformed: it should give the numbers of individuals and of loci"
    )
  }
  list(format = geno_format, n_indiv = counts[1], n_snp = counts[2])
}

# Check that a packedancestrymap genotype file has the numbers of individuals
# and of loci of its .ind and .snp files (in its header), and the size that
# they imply:
# - GENO: (n_snp + 1) records (the header and the loci) of
#   max(48, ceiling(n_indiv / 4)) bytes each
# - TGENO: a header of 48 bytes, and n_indiv records of ceiling(n_snp / 4)
#   bytes each
# A smaller file is an error (it is truncated); a larger one gives a warning
# (the extra bytes at its end are ignored). Returns the format ("GENO" or
# "TGENO").
check_packedancestry_file <- function(geno_file, n_indiv, n_snp,
                                      header = NULL) {
  if (is.null(header)) {
    header <- read_packedancestry_header(geno_file)
  }
  if (header$n_indiv != n_indiv) {
    stop(
      "Number of individuals in geno file does not match the number of ",
      "individuals in the ind file"
    )
  }
  if (header$n_snp != n_snp) {
    stop(
      "Number of variants in geno file does not match the number of ",
      "variants in the snp file"
    )
  }
  # (sizes as doubles, which are exact up to 2^53 bytes)
  if (header$format == "TGENO") {
    expected_size <- 48 + as.numeric(n_indiv) * ceiling(n_snp / 4)
  } else {
    expected_size <- (as.numeric(n_snp) + 1) * max(48, ceiling(n_indiv / 4))
  }
  geno_size <- file.size(geno_file)
  if (geno_size < expected_size) {
    stop(
      "the .geno file ", geno_file, " (", header$format, ") is smaller than ",
      "expected for the number of individuals in the .ind file (", n_indiv,
      ") and of loci in the .snp file (", n_snp, "); it might be truncated ",
      "or corrupted"
    )
  }
  if (geno_size > expected_size) {
    warning(
      "the .geno file ", geno_file, " (", header$format, ") is larger than ",
      "expected for the number of individuals in the .ind file (", n_indiv,
      ") and of loci in the .snp file (", n_snp, "); the extra bytes at its ",
      "end were ignored"
    )
  }
  header$format
}
