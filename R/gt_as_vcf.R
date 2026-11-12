#' Convert a `gen_tibble` to a VCF
#'
#' This function write a VCF from a `gen_tibble`. Phased genotypes (see
#' `phased` in [gen_tibble()]) are written as phased (e.g. `0|1`), and unphased
#' ones as unphased (e.g. `0/1`). The genotypes are read and
#' written in chunks of loci, so that this works on datasets that do not fit in
#' memory.
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
    # written in C++, one line per locus)
    gt_write_text_by_locus(
      BM = .gt_get_fbm(x),
      rowInd = .gt_fbm_rows(x),
      colInd = loci_tbl$big_index[this_chunk],
      lut = genotypes_lut,
      prefix = loci_prefix,
      sep = "\t",
      file = path.expand(file)
    )
  }
  return(file)
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
