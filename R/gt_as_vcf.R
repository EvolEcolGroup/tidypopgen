#' Convert a `gen_tibble` to a VCF
#'
#' This function write a VCF from a `gen_tibble`. Phased genotypes (see
#' `phased` in [gen_tibble()]) are written as phased (e.g. `0|1`), and unphased
#' ones as unphased (e.g. `0/1`). The genotypes are read and
#' written in chunks of loci, so that this works on datasets that do not fit in
#' memory.
#'
#' On sex chromosomes and the mitochondrial genome (see
#' [gt_chromosome_types()]), hemizygous genotypes (stored as homozygotes, as
#' in PLINK) are written as haploid genotypes (e.g. `1` rather than `1/1`):
#' those of individuals of the heterogametic sex (males with X and Y, females
#' with Z and W; from the `sex` column) on X (or Z) outside the
#' pseudo-autosomal regions and on Y (or W), and those of all individuals on
#' MT (and on Y or W if their sex is unknown). Heterozygous genotypes of
#' hemizygous individuals are written as missing (`.`), as PLINK treats them
#' as missing, with a message giving their number, and the genotypes of the
#' homogametic sex on Y (or W) are written as missing.
#'
#' @param x a [`gen_tibble`], with population coded as 'population'
#' @param file the .vcf file name with a path, or NULL (the default) to use the
#'   location of the backing files.
#' @param chunk_size the number of loci processed at a time (also used for the
#'   number of chromosomes processed at a time when writing the contig lines
#'   of the header). If left to NULL, 100,000 loci are processed at a time.
#' @param overwrite logical, should the file be overwritten if it already
#'   exists?
#' @returns the path of the .vcf file
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Write a vcf file
#' example_gt %>% gt_as_vcf()
gt_as_vcf <- function(x, file = NULL, chunk_size = NULL, overwrite = FALSE) {
  # check that x is a gen_tibble
  stopifnot_gen_tibble(x)

  # vcf writing currently only works for diploid data
  stopifnot_diploid(x$genotypes)

  if (is.null(file)) {
    file <- bigstatsr::sub_bk(
      attr(x$genotypes, "fbm")$backingfile,
      paste0(".vcf")
    )
  }
  if (file_ext(file) != "vcf") {
    file <- paste0(file, ".vcf")
  }
  if (!overwrite && file.exists(file)) {
    stop("file ", file, " already exists; use 'overwrite=TRUE' to overwrite it")
  }

  # if chunk is null, process 100,000 loci at a time (the genotypes are
  # written in C++, so the memory used by a chunk is only the information on
  # its loci)
  if (is.null(chunk_size)) {
    chunk_size <- 1e5
  }

  # set up chunks
  chunks_vec <- rep(chunk_size, floor(count_loci(x) / chunk_size))
  if (count_loci(x) %% chunk_size != 0) {
    chunks_vec <- c(chunks_vec, count_loci(x) %% chunk_size)
  }
  chunks_vec_index <- c(0, cumsum(chunks_vec))

  # generate the header (the contig lines, one per chromosome, are added
  # between its start and its end when it is written)
  vcf_header_start <- c(
    "##fileformat=VCFv4.3",
    paste0("##fileDate=", format(Sys.time(), "%Y%m%d")),
    paste0("##source=tidypopgen_v", utils::packageVersion("tidypopgen"))
  )
  vcf_header_end <- c(
    '##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">', # nolint
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste0(
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
      paste(x$id, collapse = "\t")
    )
  )
  # create copy of loci table
  loci_tbl <- show_loci(x)
  # write the header to file
  vcf_write_header(
    file = file,
    header_start = vcf_header_start,
    header_end = vcf_header_end,
    chromosome = loci_tbl$chromosome,
    position = loci_tbl$position,
    chunk_size = chunk_size
  )
  if (is_phased(x)) {
    # phased genotypes: the string of each raw value gives the alleles of the
    # two haplotypes (see phased_bytes)
    genotypes_lut <- paste0(phased_bytes[, 1], "|", phased_bytes[, 2])
    genotypes_lut[is.na(phased_bytes[, 1]) | is.na(phased_bytes[, 2])] <- ".|."
  } else {
    # create a lookup table with the string for each raw value of the FBM, by
    # recoding its decoded dosages (dosages other than 0, 1 and 2, e.g. imputed
    # dosages, are written as numbers)
    genotypes_lut <- .gt_get_fbm(x)$code256
    genotypes_lut[genotypes_lut == 0] <- "0/0"
    genotypes_lut[genotypes_lut == 1] <- "0/1"
    genotypes_lut[genotypes_lut == 2] <- "1/1"
    genotypes_lut[is.na(genotypes_lut)] <- "./."
  }
  # on sex chromosomes and mtDNA, the genotypes of hemizygous individuals
  # (e.g. males on X) are written as haploid genotypes (see vcf_haploid_luts())
  hemizygous <- vcf_haploid_luts(x, genotypes_lut)
  n_het_haploid <- 0
  # iterate over loci in chunks and append to the vcf body
  for (chunk_i in seq_along(chunks_vec)) {
    this_chunk <- (chunks_vec_index[chunk_i] + 1):chunks_vec_index[chunk_i + 1]
    # subset loci to this chunk
    loci_sub <- loci_tbl[this_chunk, ]
    # add the other columns needed for the
    loci_cols <- c("chromosome", "position", "name", "allele_ref", "allele_alt")
    loci_sub <- loci_sub %>%
      select(any_of(loci_cols)) %>%
      mutate(
        position = format_position(.data$position),
        qual = ".",
        filter = ".",
        info = "PR",
        format = "GT"
      )
    # write missing values as "." (the columns are converted to text first, as
    # a factor, e.g. of chromosomes, or integer names can not store ".")
    loci_sub <- lapply(loci_sub, function(col) {
      col <- as.character(col)
      col[is.na(col)] <- "."
      col
    })
    # the start of each line, with the information on the locus
    loci_prefix <- do.call(paste, c(unname(loci_sub), sep = "\t"))
    # append the lines of this chunk to the file (the genotypes are read and
    # written in C++, one line per locus), a run of loci of the same type at
    # a time (each type gives a set of strings to each individual)
    chunk_types <- hemizygous$locus_types[this_chunk]
    run_ends <- c(which(diff(chunk_types) != 0), length(this_chunk))
    run_starts <- c(1, run_ends[-length(run_ends)] + 1)
    for (r in seq_along(run_starts)) {
      run <- run_starts[r]:run_ends[r]
      n_het_haploid <- n_het_haploid + gt_write_text_by_locus(
        BM = .gt_get_fbm(x),
        rowInd = .gt_fbm_rows(x),
        colInd = loci_tbl$big_index[this_chunk[run]],
        lut = hemizygous$lut,
        prefix = loci_prefix[run],
        sep = "\t",
        file = path.expand(file),
        indiv_lut = hemizygous$indiv_lut[[chunk_types[run[1]]]],
        flag = hemizygous$flag
      )
    }
  }
  if (n_het_haploid > 0) {
    message(
      n_het_haploid, " heterozygous genotypes of hemizygous individuals ",
      "(e.g. males on X) were written as missing"
    )
  }
  return(file)
}

# The sets of strings used by gt_as_vcf() to write the genotypes of the
# individuals of x on each type of chromosome, from `lut`, the strings of
# diploid genotypes (256, one per raw value of the FBM): a list with `lut`
# (256 strings for each of 3 sets: diploid genotypes, haploid genotypes, and
# missing genotypes), `flag` (which strings are heterozygous haploid
# genotypes, written as missing), `locus_types` (the code of the type of each
# locus, see loci_type_index()), and `indiv_lut`, for each type, the set of
# strings of each individual (an empty vector for all diploid). Individuals
# of the heterogametic sex (males with X and Y, females with Z and W) are
# haploid on X/Z (outside the pseudo-autosomal regions) and Y/W, those of the
# other sex are diploid on X/Z and missing on Y/W, and those of unknown sex
# are diploid on X/Z and haploid on Y/W; everyone is haploid on MT. A
# homozygote (0 or 2, as PLINK stores hemizygous genotypes) is written as a
# haploid allele (0 or 1), and a heterozygote as missing (as PLINK treats
# heterozygous haploid genotypes as missing).
vcf_haploid_luts <- function(x, lut) {
  locus_types <- loci_type_index(x$genotypes)
  n_types <- length(chromosome_type_levels)
  # all individuals use the first set (diploid genotypes) by default
  indiv_lut <- rep(list(integer(0)), n_types)
  # the haploid strings, from the decoded genotypes of each raw value
  if (is_phased(x)) {
    hom <- phased_bytes[, 1] == phased_bytes[, 2]
    haploid_lut <- ifelse(!is.na(hom) & hom, phased_bytes[, 1], ".")
    het <- !is.na(hom) & !hom
  } else {
    dosage <- .gt_get_fbm(x)$code256
    haploid_lut <- lut
    haploid_lut[dosage %in% 0] <- "0"
    haploid_lut[dosage %in% 2] <- "1"
    haploid_lut[dosage %in% 1] <- "."
    haploid_lut[is.na(dosage)] <- "."
    het <- dosage %in% 1
  }
  haploid_lut <- as.character(haploid_lut)
  present <- unique(locus_types)
  hemizygous_codes <- match(hemizygous_types, chromosome_type_levels)
  if (any(present %in% hemizygous_codes)) {
    sex <- .gt_sex(x)
    sex_system <- gt_sex_system(locus_types)
    heterogametic <- if (sex_system == "ZW") 2L else 1L
    is_hetero <- sex %in% heterogametic
    is_homo <- !is.na(sex) & !is_hetero
    code_of <- function(type) match(type, chromosome_type_levels)
    # X (or Z): haploid in the heterogametic sex
    for (type in c("X", "Z")) {
      indiv_lut[[code_of(type)]] <- as.integer(is_hetero)
    }
    # Y (or W): missing in the homogametic sex, haploid otherwise
    for (type in c("Y", "W")) {
      indiv_lut[[code_of(type)]] <- ifelse(is_homo, 2L, 1L)
    }
    # MT: haploid in everyone
    indiv_lut[[code_of("MT")]] <- rep(1L, nrow(x))
  }
  list(
    lut = c(lut, haploid_lut, rep(".", 256)),
    flag = c(rep(FALSE, 256), het, rep(FALSE, 256)),
    locus_types = locus_types,
    indiv_lut = indiv_lut
  )
}


# write the header of a VCF: the lines in `header_start`, then a ##contig
# line for each chromosome, then the lines in `header_end`. Chromosomes are
# listed in the order in which they first appear (missing chromosomes last, as
# dplyr groups them, with "." as ID), with their largest position + 1 as their
# length (NA if a position of the chromosome is missing, as max() would
# give). Some data have millions of chromosomes (e.g. RAD loci assembled de
# novo), so the lengths are found without grouping or looping over
# chromosomes, and the contig lines are made and written `chunk_size`
# chromosomes at a time.
vcf_write_header <- function(
  file,
  header_start,
  header_end,
  chromosome,
  position,
  chunk_size
) {
  file_con <- file(file, open = "w")
  on.exit(close(file_con))
  writeLines(header_start, file_con)
  if (length(chromosome) > 0) {
    # the integer codes of a factor are faster to match than its labels
    chrom_key <- chromosome
    if (is.factor(chrom_key)) {
      chrom_key <- as.integer(chrom_key)
    }
    # the first locus of each chromosome, with missing chromosomes last
    first <- which(!duplicated(chrom_key))
    first <- c(first[!is.na(chrom_key[first])], first[is.na(chrom_key[first])])
    # the index of the chromosome (in the order above) of each locus
    chrom_index <- match(chrom_key, chrom_key[first])
    # sort the loci by chromosome and then by position (missing positions
    # last), so that the last locus of each chromosome has its largest
    # position
    ord <- order(chrom_index, position, method = "radix")
    is_last <- c(diff(chrom_index[ord]) != 0L, TRUE)
    chrom_max <- position[ord][is_last]
    # write the contig lines in chunks of chromosomes
    for (chunk_start in seq(1, length(first), by = chunk_size)) {
      this_chunk <- chunk_start:min(chunk_start + chunk_size - 1, length(first))
      # a missing chromosome is written as ".", as in the body of the VCF
      chrom_id <- as.character(chromosome[first[this_chunk]])
      chrom_id[is.na(chrom_id)] <- "."
      writeLines(
        paste0(
          "##contig=<ID=",
          chrom_id,
          ",length=",
          format_position(chrom_max[this_chunk] + 1),
          ">"
        ),
        file_con
      )
    }
  }
  writeLines(header_end, file_con)
}

# format positions as plain integers for text files (positions stored as
# doubles would otherwise be written in scientific notation, e.g. 1e+05)
format_position <- function(x) {
  out <- format(x, scientific = FALSE, trim = TRUE)
  out[is.na(x)] <- NA
  out
}

# TODO for tests
# read pop_a_vcf
# rewrite it as a vcf
# re-read it and check that we got back what we started with
