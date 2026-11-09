# Synthetic VCFs (GT only) to benchmark gt_merge_vcfs() and gt_concat_vcfs()
# against bcftools merge / concat followed by gen_tibble().
#   Rscript make_data.R <out_dir> [scale]  (see prepare_data.sh)
# scale multiplies the number of sites (1 = the default sizes below).
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
scale <- if (length(args) >= 2) as.numeric(args[2]) else 1
set.seed(42)
chroms <- paste0("chr", 1:22)
contig_lines <- sprintf("##contig=<ID=%s,length=250000000>", chroms)
header <- function(samples) {
  c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    contig_lines,
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", samples
    ), collapse = "\t")
  )
}
# the sites of a chromosome: sorted positions, REF and a different ALT
make_sites <- function(chr, n) {
  bases <- c("A", "C", "G", "T")
  ref <- sample.int(4, n, replace = TRUE)
  alt <- (ref + sample.int(3, n, replace = TRUE) - 1L) %% 4L + 1L
  data.table(
    CHROM = chr, POS = sort(sample.int(2.4e8, n)), ID = ".",
    REF = bases[ref], ALT = bases[alt], QUAL = ".", FILTER = "PASS",
    INFO = ".", FORMAT = "GT"
  )
}
# genotypes (dosages) of n samples at sites with alternate allele frequency p
make_dosages <- function(p, n) {
  matrix(rbinom(length(p) * n, 2, rep(p, n)), length(p), n)
}
gt_strings <- function(dosage, missing_rate = 0.02) {
  gs <- matrix(c("0/0", "0/1", "1/1")[dosage + 1L], nrow(dosage))
  gs[runif(length(gs)) < missing_rate] <- "./."
  gs
}
append_rows <- function(file, sites, gs) {
  dt <- cbind(sites, as.data.table(gs))
  fwrite(dt, file,
    sep = "\t", quote = FALSE, col.names = FALSE, append = TRUE,
    showProgress = FALSE
  )
}

# 1. concat: 100 samples, one file per chromosome (200k SNPs each), and the
#    same files with the samples in a random order in each file
n_conc <- 100
conc_samples <- sprintf("ind%03d", seq_len(n_conc))
dir.create(file.path(out_dir, "concat"), recursive = TRUE)
dir.create(file.path(out_dir, "concat_shuffled"), recursive = TRUE)
for (chr in chroms) {
  sites <- make_sites(chr, round(2e5 * scale))
  gs <- gt_strings(make_dosages(runif(nrow(sites), 0.01, 0.5), n_conc))
  f <- file.path(out_dir, "concat", paste0(chr, ".vcf"))
  writeLines(header(conc_samples), f)
  append_rows(f, sites, gs)
  ord <- sample.int(n_conc)
  f <- file.path(out_dir, "concat_shuffled", paste0(chr, ".vcf"))
  writeLines(header(conc_samples[ord]), f)
  append_rows(f, sites, gs[, ord])
}

# 2. merge, as SGDP: 50 single-sample files, each with only the sites where
#    its individual carries a non-reference allele (4M sites in all)
n_ind <- 50
ind_files <- file.path(out_dir, "merge_sgdp", sprintf("S%02d.vcf", seq_len(n_ind)))
dir.create(file.path(out_dir, "merge_sgdp"), recursive = TRUE)
for (i in seq_len(n_ind)) writeLines(header(sprintf("S%02d", i)), ind_files[i])
for (chr in chroms) {
  sites <- make_sites(chr, round(4e6 / 22 * scale))
  dosage <- make_dosages(runif(nrow(sites), 0.01, 0.5), n_ind)
  for (i in seq_len(n_ind)) {
    has_alt <- dosage[, i] > 0
    append_rows(ind_files[i], sites[has_alt],
      matrix(c("0/1", "1/1")[dosage[has_alt, i]], ncol = 1)
    )
  }
}

# 3. merge, populations: 5 files of 20 samples, genotyped at the same 2M
#    sites, with different allele frequencies
n_pop <- 5
n_per_pop <- 20
pop_files <- file.path(out_dir, "merge_pops", sprintf("pop%d.vcf", seq_len(n_pop)))
dir.create(file.path(out_dir, "merge_pops"), recursive = TRUE)
for (k in seq_len(n_pop)) {
  writeLines(header(sprintf("pop%d_%02d", k, seq_len(n_per_pop))), pop_files[k])
}
for (chr in chroms) {
  sites <- make_sites(chr, round(2e6 / 22 * scale))
  p <- runif(nrow(sites), 0.01, 0.5)
  for (k in seq_len(n_pop)) {
    p_k <- pmin(pmax(p + rnorm(length(p), 0, 0.05), 0), 1)
    append_rows(pop_files[k], sites, gt_strings(make_dosages(p_k, n_per_pop)))
  }
}
cat("done\n")
