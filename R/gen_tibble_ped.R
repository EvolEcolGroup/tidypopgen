# A function to read ped files.
# The genotypes are written straight into the backing file while the .ped file
# is read (see ped_to_fbm() in src/ped_to_fbm.cpp), so that the genotype
# matrix is never held in memory. The alleles of each locus are coded as in
# snpStats::read.pedfile() (which an earlier R reader was modified from): the
# reference allele is the first non-missing allele found in the file, and the
# alternate allele the first one that differs from it; the coding of a
# genotype is therefore known as soon as it is read, and only multi-allelic
# loci (set to missing) have to be rewritten at the end.
gen_tibble_ped <- function(
  x,
  ...,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  allow_duplicates = FALSE,
  backingfile = NULL,
  n_cores = 1,
  quiet = FALSE
) {
  # Substitute .ped with .map
  map_file <- sub_ext_ignore_case(x, "ped", ".map")
  if (!file.exists(map_file)) {
    stop("map file ", map_file, " does not exist")
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || !is.finite(n_cores) ||
        n_cores < 1) {
    stop("'n_cores' should be a positive integer")
  }

  # the number of individuals (the lines that are not blank), which the FBM
  # needs before the genotypes are written to it
  n_indiv <- count_ped_individuals(x, n_threads = n_cores)
  if (n_indiv == 0) {
    stop("Nothing read")
  }
  map <- read_ped_map(map_file)

  # create the FBM (diploid dosages, with 3 for missing genotypes), and remove
  # its files if we stop before the gen_tibble is created
  if (is.null(backingfile)) {
    backingfile <- tempfile()
  } else {
    backingfile <- change_duplicated_file_name(backingfile)
  }
  code256 <- rep(NA_real_, 256)
  code256[1:3] <- 0:2
  fbm_obj <- bigstatsr::FBM.code256(
    nrow = n_indiv,
    ncol = nrow(map),
    code = code256,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )
  fbm_files <- c(
    fbm_obj$backingfile,
    bigstatsr::sub_bk(fbm_obj$backingfile, ".rds")
  )
  fbm_done <- FALSE
  # (the FBM is dropped first, so that its backing file can be unmapped, see
  # remove_fbm_files())
  on.exit(
    if (!fbm_done) {
      rm(fbm_obj)
      remove_fbm_files(fbm_files)
    },
    add = TRUE
  )
  # sidecars next to the new backing file are left from older files with the
  # same name, and can not belong to the new data
  remove_stale_sidecars(fbm_obj$backingfile)

  # read the genotypes into the FBM. 0 is always a missing value (as in
  # check_allele_alphabet()), so that a missing genotype "0 0" (the PLINK
  # default, also written by gt_as_plink()) is never read as an allele "0"
  # when custom missing_alleles are given
  res <- ped_to_fbm(x,
    BM = fbm_obj,
    locus_names = as.character(map[, 2]),
    na_strings = union(missing_alleles, "0"),
    n_threads = n_cores
  )
  if (!quiet) {
    if (res$n_no_data > 0) {
      message("no data for ", res$n_no_data, " loci")
    }
    if (res$n_monomorphic > 0) {
      message(res$n_monomorphic, " loci were monomorphic")
    }
    if (res$n_multiallelic > 0) {
      message(res$n_multiallelic, " loci were multi-allelic --- set to NA")
    }
  }

  # the subjects need unique identifiers: the family IDs, the individual IDs,
  # or their combination (as in snpStats::read.pedfile())
  if (anyDuplicated(res$family_id) &&
        anyDuplicated(res$sample_id) &&
        anyDuplicated(paste(res$family_id, res$sample_id, sep = "."))) {
    stop("could not create unique subject identifiers")
  }
  # the first six columns of the ped file, as a fam file (IDs as characters,
  # and the type of sex and phenotype guessed, as for a fam file)
  fam <- data.frame(
    family_id = res$family_id,
    sample_id = res$sample_id,
    paternal_id = res$paternal_id,
    maternal_id = res$maternal_id,
    sex = utils::type.convert(res$sex, as.is = TRUE),
    affection = utils::type.convert(res$affection, as.is = TRUE)
  )
  # the loci table, from the map file and the alleles found in the ped file
  names(map) <- c("chromosome", "name", "genetic_dist", "position")
  loci <- map[, c("name", "chromosome", "position", "genetic_dist")]
  loci$allele_ref <- res$allele_ref
  loci$allele_alt <- res$allele_alt
  rm(map, res)

  # validate the tables (as gen_tibble.matrix() does)
  loci <- validate_loci(loci,
    check_alphabet = TRUE,
    check_duplicates = TRUE,
    allow_duplicates = allow_duplicates,
    harmonise_loci = TRUE,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles
  )
  indiv_meta <- validate_indiv_meta(indiv_meta_from_fam(fam))

  # save the FBM, and create the gen_tibble
  fbm_path <- fbm_files[2]
  saveRDS(fbm_obj, fbm_path)
  indiv_id <- indiv_meta$id
  indiv_meta <- as.list(indiv_meta)
  indiv_meta$genotypes <- new_vctrs_bigsnp(
    fbm_obj,
    fbm_file = fbm_path,
    loci = loci,
    indiv_id = indiv_id,
    ploidy = 2,
    fbm_ploidy = rep(2, n_indiv)
  )
  new_gen_tbl <- tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )
  fbm_done <- TRUE
  new_gen_tbl
}

# read the map file of a ped file (as snpStats::read.pedfile() does: the
# names of the loci are the first column of the map file without duplicates)
read_ped_map <- function(map_file) {
  map <- utils::read.table(map_file, comment.char = "")
  which <- 1
  while (anyDuplicated(map[, which])) {
    if (which == ncol(map)) {
      stop("No unambiguous snp names found on file")
    }
    which <- which + 1
  }
  map
}
