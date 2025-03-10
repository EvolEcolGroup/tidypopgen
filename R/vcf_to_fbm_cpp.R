#' Convert vcf to FBM.
#'
#' Convert a vcf file to a Filebacked Big Matrix (FBM) object. This should work
#' even for large vcf files that would not fit in memory.
#'
#' @param vcf_path the path to the vcf
#' @param chunk_size the chunk size to use on the vcf when loading the file. If
#'   NULL, a best guess will be made.
#' @param backingfile the name of the file to use as the backing file for the
#'   FBM. If NULL, the vcf path will be used.
#' @param n_cores the number of cores to use when reading the vcf file. Default
#'   is 1.
#' @return path to the resulting rds file as class bigSNP.
#' @keywords internal

vcf_to_fbm_cpp <- function(
    vcf_path,
    chunk_size = NULL,
    backingfile = NULL,
    n_cores = 1,
    quiet = FALSE) {
  if (is.null(backingfile)) {
    backingfile <- vcf_path
    backingfile <- sub("\\.vcf.gz$", "", backingfile)
    backingfile <- sub("\\.vcf$", "", backingfile)
  }
  if (file_ext(backingfile) == "bk") {
    backingfile <- bigstatsr::sub_bk(backingfile, "")
  }

  # figure out no_individuals and ploidy from the first marker
  vcf_meta <- get_ploidy_from_VCF(vcf_path)
  ploidy <- vcf_meta$ploidy
  no_individuals <- length(ploidy)
  max_ploidy <- max(ploidy)

  # if chunk is null, get the best guess of an appropriate number
  if (is.null(chunk_size)) {
    chunk_size <- bigstatsr::block_size(no_individuals)
  }

  # preassign a genotype matrix that we will then fill at each pass of the vcf
  genotypes_matrix <- matrix(NA_integer_,
    ncol = chunk_size,
    nrow = no_individuals
  )

  # set up codes for the appropriate ploidy level
  code256 <- rep(NA_real_, 256)
  code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)

  # metadata
  fam <- tibble(
    family.ID = vcf_meta$sample_names,
    sample.ID = vcf_meta$sample_names,
    paternal.ID = 0,
    maternal.ID = 0,
    sex = 0,
    affection = -9,
    ploidy = ploidy
  )

  loci <- tibble(
    chromosome = NULL,
    marker.ID = NULL,
    physical.pos = NULL,
    allele1 = NULL,
    allele2 = NULL
  )

  # create the file backed matrix
  file_backed_matrix <- bigstatsr::FBM.code256(
    nrow = no_individuals,
    ncol = 0,
    code = code256,
    backingfile = backingfile
  )

  # flag that will be changed when a chunk reaches the end of the file
  file_end <- FALSE
  start_locus <- 0

  while (!file_end) {
    # we need a function that reads chunk_size loci, staring from start_locus
    # it edits the genotype_matrix, plus returns the number of valid columns
    # (i.e. number of biallelic markers in this chunk), and a loci table
    chunk_info <- extractAltAlleleCountsFromVCF(
      vcf_path,
      genotypes_matrix, ploidy,
      no_individuals, max_ploidy + 1,
      chunk_size, start_locus
    )
    # expand the file backed matrix according to the size of the gt matrix
    # get current size
    index_start <- dim(file_backed_matrix)[2] + 1
    # add the new columns
    file_backed_matrix$add_columns(chunk_info$num_loci)
    # fill them in
    write_to_FBM(file_backed_matrix,
      allele_counts = genotypes_matrix,
      col_start = index_start - 1,
      n_loci = chunk_info$num_loci,
      ncores = n_cores
    )
    # update the metadata
    loci <- rbind(loci, chunk_info$loci_table)

    # update starting point for next chunk
    start_locus <- start_locus + chunk_size
    # check if we are at the end of file
    file_end <- chunk_info$end_of_file
  }

  # save it
  file_backed_matrix$save()

  # add an empty genetic.pos column
  loci <- loci %>% mutate(genetic.dist = 0, .after = "physical.pos")
  loci$physical.pos <- as.integer(loci$physical.pos)
  # if loci names are missing, create a name with scaffold and position
  # (same behaviour as vcfR)
  missing_loci_ids <- which(loci$marker.ID == ".")
  if (length(missing_loci_ids) > 0) {
    loci$marker.ID[missing_loci_ids] <- paste(loci$chromosome[missing_loci_ids],
      loci$physical.pos,
      sep = "_"
    )
  }

  bigsnp_obj <- structure(list(
    genotypes = file_backed_matrix,
    fam = fam,
    map = loci
  ), class = "bigSNP")

  bigsnp_obj <- bigsnpr::snp_save(bigsnp_obj)
  # and return the path to the rds
  bigsnp_obj$genotypes$rds
}

# get ploidy for a given individual
get_ploidy <- function(x) {
  max(sapply(strsplit(x, "[/|]"), function(x) length(x)))
}
