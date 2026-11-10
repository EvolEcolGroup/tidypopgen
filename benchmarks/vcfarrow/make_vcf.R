# Synthetic DiscoSnpRad VCF with the shape of the P. vaillantii dataset used
# in the VCFArrow comparison (18,579,683 variants x 56 samples, ~26 GB), built
# by resampling the 10,000 x 18 subsample that ships with VCFArrow
# (inst/extdata/vaillantii_discosnp_sub.vcf.gz).
#
#   Rscript make_vcf.R <n_variants> <out.vcf> [seed]
#
# What is kept from the real data:
# - each row copies REF, ALT, INFO and the per-sample FORMAT fields
#   (GT:DP:PL:AD:HQ) of a random row of the subsample, so the missingness
#   (~34%), the 7% insertions and the line length (~1.5 kB with 56 samples)
#   follow the real data;
# - each of the 56 samples takes, at each row, the fields of a random sample
#   of the same lineage (GS, BS, WA, outgroup OG), so allele frequencies
#   differ between groups;
# - CHROM is a DiscoSnpRad path ("SNP_higher_path_<n>"), and paths carry 1-5
#   SNPs in the proportions implied by the ID suffixes of the subsample, so
#   ~74% of rows start a new "chromosome" (~13.8 M chromosomes at full size);
# - rows are in descending lexicographic order of CHROM, as in the subsample.
# Positions within a path are random (30-160), strictly increasing within the
# path.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
n_variants <- as.numeric(args[1])
out_file <- args[2]
seed <- if (length(args) >= 3) as.integer(args[3]) else 1L
chunk_rows <- 250000L

here <- dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)
)))
src_vcf <- file.path(here, "data", "vaillantii_discosnp_sub.vcf.gz")
src_strata <- file.path(here, "data", "strata")

set.seed(seed)

# source data ---------------------------------------------------------------
hdr <- readLines(src_vcf, n = 200)
meta_lines <- grep("^##", hdr, value = TRUE)
src <- fread(src_vcf, skip = "#CHROM", sep = "\t", header = TRUE,
  colClasses = "character", showProgress = FALSE)
src_samples <- names(src)[-(1:9)]
src_fields <- as.matrix(src[, src_samples, with = FALSE])
strata <- fread(src_strata)
src_group <- strata$pop[match(src_samples, strata$id)]
stopifnot(!anyNA(src_group))

# 56 output samples, by lineage
group_sizes <- c(GS = 18L, BS = 18L, WA = 17L, OG = 3L)
out_group <- rep(names(group_sizes), group_sizes)
out_samples <- unlist(lapply(names(group_sizes), function(g) {
  sprintf("%s_%02d", g, seq_len(group_sizes[[g]]))
}))
# columns of the source matrix each output sample draws from
src_cols_by_group <- split(seq_along(src_samples), src_group)

# path structure --------------------------------------------------------------
# paths with 1, 2, 3, 4 and 5 SNPs, from the ID suffixes of the subsample
path_size_prob <- c(5878, 918, 372, 170, 98)
mean_size <- sum(seq_along(path_size_prob) * path_size_prob) /
  sum(path_size_prob)
n_paths <- ceiling(n_variants / mean_size * 1.02)
path_size <- sample.int(5L, n_paths, replace = TRUE, prob = path_size_prob)
n_paths <- which(cumsum(path_size) >= n_variants)[1]
path_size <- path_size[seq_len(n_paths)]
path_size[n_paths] <- path_size[n_paths] - (sum(path_size) - n_variants)

path_id <- sample.int(max(2e7, ceiling(1.8 * n_paths)), n_paths)
# descending lexicographic order of the decimal strings, without making them:
# left-align the digits (key); among equal keys the longer string is larger
n_digits <- floor(log10(path_id)) + 1
key <- path_id * 10^(9 - n_digits)
ord <- order(key, n_digits, decreasing = TRUE, method = "radix")
path_id <- path_id[ord]
path_size <- path_size[ord]
rm(key, n_digits, ord)

# write ------------------------------------------------------------------------
meta_lines <- c(
  meta_lines[1],
  paste0("##source=tidypopgen_perf/vcfarrow/make_vcf.R (synthetic, seed ",
    seed, ", resampled from VCFArrow vaillantii_discosnp_sub.vcf.gz)"),
  meta_lines[-c(1, grep("^##source=", meta_lines))]
)
writeLines(c(meta_lines, paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
  "FILTER", "INFO", "FORMAT", out_samples), collapse = "\t")), out_file)

# chunks of ~chunk_rows rows, ending at path boundaries
path_chunk <- (cumsum(path_size) - 1) %/% chunk_rows
t0 <- Sys.time()
for (p in split(seq_len(n_paths), path_chunk)) {
  sizes <- path_size[p]
  n <- sum(sizes)
  k <- sequence(sizes) # SNP index within path
  pid <- rep(path_id[p], sizes)
  id <- ifelse(rep(sizes, sizes) == 1L, as.character(pid),
    paste0(pid, "_", k))
  # strictly increasing positions within path: a start (30-100) plus gaps
  # of 1-15 bp (no duplicated chromosome-position pairs)
  gap <- sample.int(15L, n, replace = TRUE)
  gap[k == 1L] <- 0L
  cs <- cumsum(gap)
  pos <- rep(sample(30:100, length(p), replace = TRUE), sizes) + cs -
    rep(cs[k == 1L], sizes)
  r <- sample.int(nrow(src), n, replace = TRUE)
  fixed <- data.table(
    CHROM = paste0("SNP_higher_path_", pid), POS = pos, ID = id,
    REF = src$REF[r], ALT = src$ALT[r], QUAL = ".", FILTER = ".",
    INFO = src$INFO[r], FORMAT = src$FORMAT[r]
  )
  geno <- lapply(seq_along(out_samples), function(j) {
    cols <- src_cols_by_group[[out_group[j]]]
    src_fields[cbind(r, cols[sample.int(length(cols), n, replace = TRUE)])]
  })
  setDT(geno)
  fwrite(cbind(fixed, geno), out_file, append = TRUE, sep = "\t",
    quote = FALSE, col.names = FALSE, nThread = 2)
}
message(sprintf("wrote %s variants x %d samples to %s in %.1f s",
  format(n_variants, big.mark = ","), length(out_samples), out_file,
  as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# strata file next to the VCF (id, pop), for the "set groups" step
fwrite(data.table(id = out_samples, pop = out_group),
  sub("\\.vcf(\\.gz)?$", ".strata.tsv", out_file), sep = "\t")

# number of variants, read by bench.py (counting the lines takes minutes)
writeLines(format(n_variants, scientific = FALSE),
  sub("\\.vcf(\\.gz)?$", ".n_variants", out_file))
