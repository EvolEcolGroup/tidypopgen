# Benchmark the cpp VCF parser of two installed versions of tidypopgen
#
# Install each version in its own library, e.g. to compare the current
# branch with dev (from the root of the repository):
#   git worktree add ../tpg_dev origin/dev
#   mkdir -p ../lib_dev ../lib_new
#   R CMD INSTALL --no-docs -l ../lib_dev ../tpg_dev
#   R CMD INSTALL --no-docs -l ../lib_new .
#
# Create a test vcf with n_variants and n_samples. ploidy is cycled across
# samples (e.g. "2", "4" or "2,4"); if missing_sample is 1, the first sample
# is a bare "." on every variant (the slowest case for estimating ploidy):
#   Rscript benchmark_vcf_parser.R make big.vcf 10000000 20 2 0
#
# Time each version in its own R session, in a balanced order (a version
# that always runs in the same position can look faster or slower than it
# is):
#   for v in dev new new dev dev new; do
#     Rscript benchmark_vcf_parser.R run ../lib_$v big.vcf $v
#   done
#
# Each run prints: label, file, seconds for vcf_loci_table(), seconds for
# vcf_genotypes_to_fbm(), number of loci, number of biallelic loci, ploidies
# (ploidy x number of individuals) and the md5 of the backing file. The md5
# must be the same for both versions, unless the change is meant to alter the
# genotypes. Do not run long benchmarks on a laptop that may suspend, as the
# elapsed times include the time spent suspended.

args <- commandArgs(trailingOnly = TRUE)

# write a vcf with random genotypes (5% fully missing, 30% phased), and 2%
# multiallelic variants (which the parser does not store)
make_vcf <- function(vcf_path, n_variants, n_samples, ploidy = "2",
                     missing_sample = FALSE, chunk_size = 1e5) {
  set.seed(42)
  ploidy <- rep_len(as.integer(strsplit(ploidy, ",")[[1]]), n_samples)
  con <- file(vcf_path, open = "w")
  on.exit(close(con))
  writeLines(c(
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", paste0("ind", seq_len(n_samples))
    ), collapse = "\t")
  ), con)
  bases <- c("A", "C", "G", "T")
  for (start in seq(1, n_variants, by = chunk_size)) {
    m <- min(chunk_size, n_variants - start + 1)
    locus <- as.integer(start + seq_len(m) - 1)
    alt_freq <- stats::runif(m, 0.05, 0.95)
    # one column of genotypes (with depth) per sample
    geno <- matrix(vapply(seq_len(n_samples), function(i) {
      sep <- ifelse(stats::runif(m) < 0.3, "|", "/")
      gt <- as.character(as.integer(stats::runif(m) < alt_freq))
      for (a in seq_len(ploidy[i] - 1)) {
        gt <- paste0(gt, sep, as.integer(stats::runif(m) < alt_freq))
      }
      gt[stats::runif(m) < 0.05] <- paste(rep(".", ploidy[i]), collapse = "/")
      paste0(gt, ":", sample.int(60, m, replace = TRUE))
    }, character(m)), nrow = m)
    if (missing_sample) {
      # a bare "." gives no information on ploidy
      geno[, 1] <- "."
    }
    ref <- sample.int(4, m, replace = TRUE)
    alt <- (ref - 1 + sample.int(3, m, replace = TRUE)) %% 4 + 1
    alt_str <- bases[alt]
    multi <- which(locus %% 50 == 7)
    for (j in multi) {
      other <- setdiff(1:4, c(ref[j], alt[j]))[1]
      alt_str[j] <- paste0(alt_str[j], ",", bases[other])
    }
    chromosome <- as.integer(1 + ((locus - 1) * 22) %/% n_variants)
    meta <- paste(chromosome, locus * 100L, paste0("snp", locus), bases[ref],
      alt_str, ".", "PASS", ".", "GT:DP",
      sep = "\t"
    )
    writeLines(do.call(paste, c(list(meta), as.data.frame(geno),
      sep = "\t"
    )), con)
  }
}

# time the two steps of the cpp parser (as used by gen_tibble())
run_parser <- function(lib, vcf_path, label) {
  suppressPackageStartupMessages(library(tidypopgen, lib.loc = lib))
  time_loci <- system.time(
    vcf_meta <- tidypopgen:::vcf_loci_table(vcf_path)
  )[["elapsed"]]
  max_ploidy <- max(vcf_meta$ploidy)
  code256 <- rep(NA_real_, 256)
  code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)
  backingfile <- tempfile()
  fbm <- bigstatsr::FBM.code256(
    nrow = length(vcf_meta$ploidy),
    ncol = sum(vcf_meta$biallelic),
    code = code256,
    backingfile = backingfile
  )
  time_geno <- system.time(
    tidypopgen:::vcf_genotypes_to_fbm(vcf_path, fbm,
      biallelic = vcf_meta$biallelic,
      missing_value = max_ploidy + 1,
      n_header_lines = vcf_meta$n_header_lines
    )
  )[["elapsed"]]
  ploidy_counts <- table(vcf_meta$ploidy)
  cat(sprintf(
    "%s\t%s\t%.2f\t%.2f\t%d\t%d\t%s\t%s\n", label, basename(vcf_path),
    time_loci, time_geno, nrow(vcf_meta$loci_tbl), sum(vcf_meta$biallelic),
    paste(names(ploidy_counts), ploidy_counts, sep = "x", collapse = ","),
    unname(tools::md5sum(paste0(backingfile, ".bk")))
  ))
  unlink(paste0(backingfile, c(".bk", ".rds")))
}

if (length(args) >= 4 && args[1] == "make") {
  make_vcf(args[2], as.numeric(args[3]), as.integer(args[4]),
    ploidy = if (length(args) >= 5) args[5] else "2",
    missing_sample = length(args) >= 6 && args[6] == "1"
  )
} else if (length(args) == 4 && args[1] == "run") {
  run_parser(args[2], args[3], args[4])
} else {
  stop("see the comments at the top of this script for its usage")
}
