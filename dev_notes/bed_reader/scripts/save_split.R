# How long saving and loading a loci table of 91M rows takes, by format
set.seed(20261120)
n <- 91162139
out <- file.path(commandArgs(trailingOnly = TRUE)[1], "loci")
nuc <- c("A", "C", "G", "T")
loci <- tibble::tibble(
  big_index = seq_len(n), name = seq_len(n),
  chromosome = factor(sort(sample.int(20, n, replace = TRUE))),
  position = as.integer(cumsum(sample.int(200, n, replace = TRUE))),
  genetic_dist = 0,
  allele_ref = nuc[sample.int(4, n, replace = TRUE)],
  allele_alt = nuc[sample.int(4, n, replace = TRUE)]
)
tm <- function(label, expr) {
  t0 <- proc.time(); res <- force(expr); t1 <- proc.time() - t0
  cat(sprintf("%-40s elapsed %6.2f  size %s\n", label, t1[["elapsed"]],
              if (file.exists(paste0(out, ".rds")))
                format(file.size(paste0(out, ".rds")) / 1e9, digits = 3) else ""))
  invisible(res)
}
save_gz1 <- function(x) {
  con <- gzfile(paste0(out, ".rds"), open = "wb", compression = 1)
  on.exit(close(con)); saveRDS(x, con)
}
tm("save gzip 1 (current)", save_gz1(loci))
tm("load gzip 1 (current)", y <- readRDS(paste0(out, ".rds")))
tm("serialize to memory only (no file)", z <- serialize(loci, NULL)); rm(z)
tm("save uncompressed", saveRDS(loci, paste0(out, ".rds"), compress = FALSE))
tm("load uncompressed", y <- readRDS(paste0(out, ".rds")))
# alleles as factors (integer codes)
loci_f <- loci
loci_f$allele_ref <- factor(loci_f$allele_ref)
loci_f$allele_alt <- factor(loci_f$allele_alt)
tm("save gzip 1, factor alleles", save_gz1(loci_f))
tm("load gzip 1, factor alleles", y <- readRDS(paste0(out, ".rds")))
tm("  then as.character() the alleles",
   {y$allele_ref <- as.character(y$allele_ref); y$allele_alt <- as.character(y$allele_alt)})
unlink(paste0(out, ".rds"))
