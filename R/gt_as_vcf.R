#' Convert a `gen_tibble` to a VCF
#'
#' This function write a VCF from a `gen_tibble`. Phased genotypes (see
#' [gen_tibble_phased()]) are written as phased (e.g. `0|1`), and unphased
#' ones as unphased (e.g. `0/1`). The genotypes are read and
#' written in chunks of loci, so that this works on datasets that do not fit in
#' memory.
#'
#' @param x a [`gen_tibble`], with population coded as 'population'
#' @param file the .vcf file name with a path, or NULL (the default) to use the
#'   location of the backing files.
#' @param chunk_size the number of loci processed at a time. If left to NULL,
#'   100,000 loci are processed at a time.
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

  # generate the header
  vcf_header <- c(
    "##fileformat=VCFv4.3",
    paste0("##fileDate=", format(Sys.time(), "%Y%m%d")),
    paste0("##source=tidypopgen_v", utils::packageVersion("tidypopgen"))
  )
  # create copy of loci table
  loci_tbl <- show_loci(x)
  # reorder chromosomes levels in the order in which they appear
  loci_tbl$chromosome <- fct_inorder_base(loci_tbl$chromosome)
  # get max position for each chromosome
  chromosome_summary <- loci_tbl %>%
    group_by(.data$chromosome) %>%
    summarise(max = max(.data$position))
  for (chrom_i in seq_len(nrow(chromosome_summary))) {
    vcf_header <- c(
      vcf_header,
      paste0(
        "##contig=<ID=",
        chromosome_summary$chromosome[chrom_i],
        ",length=",
        format_position(chromosome_summary$max[chrom_i] + 1),
        ">"
      )
    )
  }
  vcf_header <- c(
    vcf_header,
    '##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">' # nolint
  )
  vcf_header <- c(
    vcf_header,
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  )
  vcf_header <- c(
    vcf_header,
    paste0(
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
      paste(x$id, collapse = "\t")
    )
  )
  # write it to file
  utils::write.table(
    vcf_header,
    file = file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
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
    loci_sub[is.na(loci_sub)] <- "."
    # the start of each line, with the information on the locus (the factor
    # of chromosomes is pasted as its labels)
    loci_prefix <- do.call(paste, c(unname(as.list(loci_sub)), sep = "\t"))
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
