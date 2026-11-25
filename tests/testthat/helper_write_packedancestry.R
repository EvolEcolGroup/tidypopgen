# Write a packedancestrymap genotype file (as EIGENSOFT's geno.c does), from
# a matrix of packed values (individuals in rows, loci in columns): the number
# of copies of the first allele of the .snp file (0, 1 or 2), or 3 if missing.
# With transposed = FALSE, a GENO file (a header record, then a record per
# locus, of max(48, ceiling(n / 4)) bytes each); with transposed = TRUE, a
# TGENO file (a header of 48 bytes, then a record per individual, of
# ceiling(m / 4) bytes each). Each byte holds 4 values, the first one in the
# highest bits; the hashes in the header are 0.
write_packed_geno <- function(values, file, transposed = FALSE) {
  n <- nrow(values)
  m <- ncol(values)
  if (transposed) {
    records <- t(values)
    record_bytes <- ceiling(m / 4)
    header_bytes <- 48
    magic <- "TGENO"
  } else {
    records <- values
    record_bytes <- max(48, ceiling(n / 4))
    header_bytes <- record_bytes
    magic <- "GENO"
  }
  # one record per column, padded with zeros
  padded <- matrix(0L, nrow = 4 * record_bytes, ncol = ncol(records))
  padded[seq_len(nrow(records)), ] <- records
  quads <- array(padded, dim = c(4, record_bytes, ncol(records)))
  bytes <- quads[1, , ] * 64L + quads[2, , ] * 16L + quads[3, , ] * 4L +
    quads[4, , ]
  header <- charToRaw(sprintf("%s %7d %7d %x %x", magic, n, m, 0L, 0L))
  header <- c(header, raw(header_bytes - length(header)))
  writeBin(c(header, as.raw(bytes)), file)
  invisible(file)
}

# Write the .snp and .ind files that go with a packedancestrymap genotype file
# of n individuals and m loci (loci on chromosome 1, alleles A and G)
write_packed_snp_ind <- function(prefix, n, m) {
  snp <- data.frame(
    name = paste0("rs", seq_len(m)),
    chromosome = 1,
    genetic_dist = 0,
    position = seq_len(m) * 100,
    allele_1 = "A",
    allele_2 = "G"
  )
  utils::write.table(snp, paste0(prefix, ".snp"),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  ind <- data.frame(
    id = paste0("ind", seq_len(n)),
    sex = "U",
    population = rep(c("pop1", "pop2"), length.out = n)
  )
  utils::write.table(ind, paste0(prefix, ".ind"),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  invisible(prefix)
}
