
gen_tibble_bed_rds <- function(
    x,
    ...,
    valid_alleles = c("A", "T", "C", "G"),
    missing_alleles = c("0", "."),
    backingfile = NULL,
    quiet = FALSE) {
  # if it is a bed file, we convert it to a bigsnpr
  if (tolower(file_ext(x)) == "bed") {
    if (is.null(backingfile)) {
      backingfile <- bigsnpr::sub_bed(x)
    }
    bigsnp_path <- bigsnpr::snp_readBed(x, backingfile = backingfile)
  } else if (tolower(file_ext(x)) == "rds") {
    bigsnp_path <- x
  }
  
  bigsnp_obj <- bigsnpr::snp_attach(bigsnp_path)
  
  if (all(bigsnp_obj$fam$family.ID == bigsnp_obj$fam$sample.ID)) {
    indiv_meta <- list(id = bigsnp_obj$fam$sample.ID)
  } else {
    indiv_meta <- list(
      id = bigsnp_obj$fam$sample.ID,
      population = bigsnp_obj$fam$family.ID
    )
  }
  # check if the bignsp_obj$fam table has ploidy column, if not, set ploidy to 2
  if ("ploidy" %in% names(bigsnp_obj$fam)) {
    ploidy <- bigsnp_obj$fam$ploidy
  } else {
    ploidy <- 2
    bigsnp_obj$fam$ploidy <- 2
  }
  
  indiv_meta$genotypes <- new_vctrs_bigsnp(
    bigsnp_obj,
    bigsnp_file = bigsnp_path,
    indiv_id = bigsnp_obj$fam$sample.ID,
    ploidy = ploidy
  )
  
  # transfer some of the fam info to the metadata table if it is not missing
  # (0 is the default missing value)
  fam_info <- .gt_get_bigsnp(indiv_meta)$fam
  if (!all(fam_info$paternal.ID == 0)) {
    indiv_meta$paternal_ID <- fam_info$paternal.ID
    indiv_meta$paternal_ID[indiv_meta$paternal_ID == 0] <- NA
  }
  if (!all(fam_info$maternal.ID == 0)) {
    indiv_meta$maternal_ID <- fam_info$maternal.ID
    indiv_meta$maternal_ID[indiv_meta$maternal_ID == 0] <- NA
  }
  # if sex is numeric
  if (inherits(fam_info$sex, "numeric")) {
    if (!all(fam_info$sex == 0)) {
      indiv_meta$sex <- dplyr::case_match(
        fam_info$sex,
        1 ~ "male",
        2 ~ "female",
        .default = NA,
        .ptype = factor(levels = c("female", "male"))
      )
    }
  }
  
  if (inherits(fam_info$affection, "numeric")) {
    if (!all(fam_info$affection %in% c(0, -9))) {
      indiv_meta$phenotype <- dplyr::case_match(
        fam_info$affection,
        1 ~ "control",
        2 ~ "case",
        -9 ~ NA,
        .default = NA,
        .ptype = factor(levels = c("control", "case"))
      )
    }
  }
  
  new_gen_tbl <- tibble::new_tibble(
    indiv_meta,
    class = "gen_tbl"
  )
  check_allele_alphabet(
    new_gen_tbl,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles,
    remove_on_fail = TRUE
  )
  show_loci(new_gen_tbl) <- harmonise_missing_values(
    show_loci(new_gen_tbl),
    missing_alleles = missing_alleles
  )
  return(new_gen_tbl)
}




#' Generate a tibble of loci from a PLINK .bim file
#'
#' This function reads a PLINK .bim file and extracts the relevant information
#' to create a loci tibble compatible with a `gen_tibble` object.
#' @param bim the path to a bim file or a data.frame containing the contents of a
#'  bim file.
#' @returns A tibble with columns: `big_index`, `name`, `chromosome`,
#' `position`, `gentic_dist`, `allele_ref`, `allele_alt`
#' @keywords internal
#' @noRd

loci_from_bim <- function(bim) {
  if (is.character(bim)) {
    bim <- read.table(bim, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(bim) || ncol(bim) < 6) {
    stop("bim must be a data.frame with at least 6 columns or a path to a .bim file")
  }
  loci <- tibble::tibble(
    big_index = seq_len(nrow(bim)),
    name = as.character(bim[[2]]),
    chromosome = as.factor(bim[[1]]),
    position = as.integer(bim[[4]]),
    genetic_dist = as.numeric(bim[[3]]),
    allele_ref = as.character(bim[[5]]),
    allele_alt = as.character(bim[[6]])
  )
  return(loci)
}
