# write a small VCF (with GT only) for the tests: loci is a data.frame with
# chrom, pos, id, ref and alt, and genotypes a character matrix (one row per
# locus, one column per sample)
write_merge_vcf <- function(samples, loci, genotypes,
                            file = tempfile(fileext = ".vcf")) {
  header <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", samples
    ), collapse = "\t")
  )
  body <- paste(loci$chrom, loci$pos, loci$id, loci$ref, loci$alt, ".",
    "PASS", ".", "GT", apply(genotypes, 1, paste, collapse = "\t"),
    sep = "\t"
  )
  writeLines(c(header, body), file)
  file
}
