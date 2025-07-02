#' Export a `gen_tibble` object to PLINK bed format
#'
#' This function exports all the information of a `gen_tibble` object into a
#' PLINK bed, ped or raw file (and associated files, i.e. .bim and .fam for
#' .bed; .fam for .ped).
#'
#' If the gen_tibble has been read in from vcf format, family.ID in the
#' resulting plink files will be the same as sample.ID. If the gen_tibble has a
#' grouping variable, this will be used as the family.ID in the resulting plink
#' files. NOTE that writing to bed has been optimised for speed, but writing to
#' ped or raw is slower, especially for large datasets.
#'
#' @param x a [`gen_tibble`] object
#' @param file a character string giving the path to output file. If left to
#'   NULL, the output file will have the same path and prefix of the
#'   backingfile.
#' @param type one of "bed", "ped" or "raw"
#' @param overwrite boolean whether to overwrite the file.
#' @param chromosomes_as_int boolean whether to use the integer representation
#'   of the chromosomes
#' @returns the path of the saved file
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Write a bed file
#' example_gt %>% gt_as_plink(type = "bed", file = paste0(tempfile(), "_plink"))
#'
#' # Write a ped file
#' example_gt %>% gt_as_plink(type = "ped", file = paste0(tempfile(), "_plink"))
#'
#' # Write a raw file
#' example_gt %>% gt_as_plink(type = "raw", file = paste0(tempfile(), "_plink"))
gt_as_plink <- function(
    x,
    file = NULL,
    type = c("bed", "ped", "raw"),
    overwrite = TRUE,
    chromosomes_as_int = FALSE) {
  # check that x is a gen_tibble
  if (!methods::is(x, "gen_tbl")) {
    stop("x must be a gen_tibble")
  }

  if (is.null(chromosomes_as_int)) {
    chromosomes_as_int <- FALSE
  }

  type <- match.arg(type)

  if (is.null(file)) {
    file <- bigstatsr::sub_bk(
      attr(
        x$genotypes,
        "bigsnp"
      )$genotypes$backingfile,
      paste0(".", type)
    )
  }
  if (file_ext(file) != type) {
    file <- paste0(file, ".", type)
  }

  if (type == "bed") {
    all_files <- c(
      file,
      bigsnpr::sub_bed(file, ".bim"),
      bigsnpr::sub_bed(file, ".fam")
    )
  } else if (type == "ped") {
    all_files <- c(
      file,
      gsub(".ped", ".map", file)
    )
  } else if (type == "raw") {
    all_files <- file
  }

  if (any(file.exists(all_files))) {
    if (overwrite) {
      file.remove(all_files[file.exists(all_files)])
    } else {
      stop(
        "at least one of",
        all_files,
        " already exists; remove if first or set 'overwrite' = TRUE"
      )
    }
  }

  if (type == "bed") {
    gt_write_bed(x, file, chromosomes_as_int)
  } else {
    gt_write_ped_raw(x, file, type)
  }
}


gt_write_bed <- function(x, file, chromosomes_as_int) {
  bed_path <- bigsnpr::snp_writeBed(
    attr(x$genotypes, "bigsnp"),
    bedfile = file,
    ind.row = vctrs::vec_data(x$genotypes),
    ind.col = show_loci(x)$big_index
  )
  # the bim and fam files written by bigsnpr contain the information of the
  # original bigsnpr object. We now update that information with the info
  # from the gen_tibble
  if (!chromosomes_as_int) {
    bim_path <- bigsnpr::sub_bed(bed_path, ".bim")
    bim_table <- show_loci(x) %>%
      dplyr::select(dplyr::all_of(c(
        "chromosome",
        "name",
        "genetic_dist",
        "position",
        "allele_alt",
        "allele_ref"
      )))
    colnames(bim_table) <- c(
      "chromosome",
      "name",
      "genetic_dist",
      "position",
      "allele_ref",
      "allele_alt"
    )
    bim_table$allele_alt[is.na(bim_table$allele_alt)] <- "0"
    bim_table$allele_ref[is.na(bim_table$allele_ref)] <- "0"
    utils::write.table(
      bim_table,
      bim_path,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  } else {
    bim_path <- bigsnpr::sub_bed(bed_path, ".bim")
    bim_table <- show_loci(x) %>%
      dplyr::select(dplyr::all_of(c(
        "chr_int",
        "name",
        "genetic_dist",
        "position",
        "allele_alt",
        "allele_ref"
      )))
    colnames(bim_table) <- c(
      "chr_int",
      "name",
      "genetic_dist",
      "position",
      "allele_ref",
      "allele_alt"
    )
    bim_table$allele_alt[is.na(bim_table$allele_alt)] <- "0"
    bim_table$allele_ref[is.na(bim_table$allele_ref)] <- "0"
    utils::write.table(
      bim_table,
      bim_path,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
  fam_path <- bigsnpr::sub_bed(bed_path, ".fam")
  fam_table <- utils::read.table(fam_path)
  fam_table$V2 <- x$id
  # if this is a group tibble, use the grouping variable
  if (inherits(x, "grouped_gen_tbl")) {
    fam_table$V1 <- x %>%
      select(dplyr::group_vars(x)) %>%
      dplyr::pull(1)
  }
  utils::write.table(
    fam_table,
    fam_path,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  # return the path to the file
  return(bed_path) # nolint
}


gt_write_ped_raw <- function(
    x,
    file,
    plink_format = c("raw", "ped"),
    chunk_size = 10000) {
  # loci information
  loci <- show_loci(x)
  # replace missing value with zero
  loci$allele_alt[is.na(loci$allele_alt)] <- "0"

  # create col names to use in the raw file
  raw_col_names <- c(
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PHENOTYPE",
    paste0(
      loci$name,
      "_",
      toupper(loci$allele_alt),
      "(/",
      toupper(loci$allele_ref),
      ")"
    )
  )

  # create the info for the fam file
  indiv_meta <- tibble(
    population = x$population,
    id = x$id,
    pat = pull_na(x, "pat"),
    mat = pull_na(x, "mat"),
    sex = pull_na(x, "sex"),
    phenotype = pull_na(x, "phenotype")
  )
  # recode some variables
  indiv_meta$sex <- dplyr::case_match(
    as.character(indiv_meta$sex),
    "female" ~ "2",
    "male" ~ "1",
    .default = "0"
  )
  indiv_meta$pat[is.na(indiv_meta$pat)] <- 0
  indiv_meta$mat[is.na(indiv_meta$mat)] <- 0
  indiv_meta$phenotype <- dplyr::case_match(
    as.character(indiv_meta$phenotype),
    "control" ~ "1",
    "case" ~ "2",
    .default = indiv_meta$phenotype
  )
  indiv_meta$phenotype[is.na(indiv_meta$phenotype)] <- -9

  # create chunks
  n_ind <- nrow(x)
  chunks <- split(1:n_ind, ceiling(seq_along(1:n_ind) / chunk_size))

  # loop over to write
  for (i_chunk in seq_len(length(chunks))) {
    chunk <- chunks[[i_chunk]]
    raw_table <- cbind(
      indiv_meta[chunk, ],
      show_genotypes(x[chunk, ])
    )

    # now recode the genotypes with letters if raw
    if (plink_format == "ped") {
      # TODO this would be faster as an apply function
      for (i in 1:(ncol(raw_table) - 6)) {
        raw_table[, i + 6] <- recode_genotype(
          raw_table[, i + 6],
          loci$allele_ref[i],
          loci$allele_alt[i]
        )
      }
    }

    colnames(raw_table) <- raw_col_names
    # append column names only the first time, when the file does not exist
    utils::write.table(
      raw_table,
      file = file,
      sep = " ",
      row.names = FALSE,
      col.names = (!file.exists(file) & plink_format == "raw"),
      append = file.exists(file),
      quote = FALSE
    )
  }

  ## create info for the map file
  loci_meta <- show_loci(x)
  map_file <- paste0(substr(file, 1, nchar(file) - 4), ".map")
  map_table <- tibble(
    chrom = pull_na(loci_meta, "chromosome"),
    id = pull_na(loci_meta, "name"),
    cM = pull_na(loci_meta, "cM"),
    position = pull_na(loci_meta, "position")
  )
  map_table[is.na(map_table)] <- 0

  utils::write.table(
    map_table,
    file = map_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  return(file)
}


# internal function returning either a vector from a given column
# or NA if it does not exist
pull_na <- function(.x, .col_name) {
  if (.col_name %in% names(.x)) {
    return(.x %>% pull(.col_name)) # nolint
  } else {
    return(rep(NA, nrow(.x))) # nolint
  }
}

# recode dosage as letter genotypes
recode_genotype <- function(x, allele_ref, allele_alt) {
  x <- as.character(x)
  genotypes <- c( # nolint
    paste(allele_ref, allele_ref),
    paste(allele_ref, allele_alt),
    paste(allele_alt, allele_alt)
  )
  x <- dplyr::case_match(
    x,
    "0" ~ genotypes[1],
    "1" ~ genotypes[2],
    "2" ~ genotypes[3]
  )
  x[is.na(x)] <- c("0 0")
  x
}
