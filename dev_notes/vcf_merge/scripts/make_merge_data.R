# Synthetic VCFs to profile gt_merge_vcfs() (files with different samples).
#   Rscript make_merge_data.R <out_dir> <scenario>
# scenarios:
#   pops : 4 files of 100 samples (GT only), ~1.05M candidate sites; each file
#          has ~97% of the sites; 1% of the loci of files 2-4 have swapped
#          REF/ALT (as in VCFs exported by PLINK), 0.1% a different ALT
#          (merge with ref_mismatch = "swap")
#   sgdp : 50 files of 1 sample, 2.5M candidate sites, each file has the
#          sites where its sample carries an ALT allele (~45%); lines carry
#          an INFO field (~150 bytes) and FORMAT GT:DP:GQ, as the per-sample
#          VCFs of a variant caller (merge with variants = "union",
#          missing_to_ref = TRUE)
#   many : 200 files of 1 sample, GT only, 1M candidate sites (~45% per
#          file); isolates the cost per file and of scattered writes
# The VCFs are written uncompressed (compress them with bgzip afterwards).
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
scenario <- args[2]
set.seed(20261119)
chroms <- as.character(1:22)
contig_lines <- sprintf("##contig=<ID=%s,length=250000000>", chroms)
header <- function(samples, info = FALSE) {
  c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    if (info) {
      c(
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">",
        # (bcftools needs the INFO fields in the header)
        sprintf("##INFO=<ID=%s,Number=1,Type=%s,Description=\"%s\">",
          c("DP", "MQ", "FS", "QD", "CpG", "RepeatMasker", "Map20",
            "TandemRepeat", "SegDup", "ANN"),
          c("Integer", "Float", "Float", "Float", "Integer", "String",
            "Integer", "Integer", "Integer", "String"),
          "annotation"
        )
      )
    },
    contig_lines,
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", samples
    ), collapse = "\t")
  )
}
# the sites of a chromosome: sorted positions, REF and a different ALT, and
# (optionally) an INFO field of per-site annotations
make_sites <- function(chr, n, info = FALSE) {
  bases <- c("A", "C", "G", "T")
  ref <- sample.int(4, n, replace = TRUE)
  alt <- (ref + sample.int(3, n, replace = TRUE) - 1L) %% 4L + 1L
  pos <- sort(sample.int(2.4e8, n))
  dt <- data.table(
    CHROM = chr, POS = pos,
    ID = ifelse(runif(n) < 0.6, paste0("rs", sample.int(9e8, n)), "."),
    REF = bases[ref], ALT = bases[alt], QUAL = ".", FILTER = "PASS",
    INFO = ".", FORMAT = "GT"
  )
  if (info) {
    dt[, INFO := sprintf(
      paste0(
        "DP=%d;MQ=%.1f;FS=%.3f;QD=%.2f;CpG=%d;RepeatMasker=%s;Map20=%d;",
        "TandemRepeat=%d;SegDup=%d;ANN=%s|intergenic_region|MODIFIER|%s"
      ),
      rpois(n, 30), runif(n, 40, 60), rexp(n), runif(n, 0, 40),
      rbinom(n, 1, 0.1), sample(c("L1", "Alu", "MIR", "."), n, TRUE),
      rbinom(n, 1, 0.9), rbinom(n, 1, 0.05), rbinom(n, 1, 0.02), ALT,
      sprintf("ENSG%011d", sample.int(6e4, n, TRUE))
    )]
    dt[, FORMAT := "GT:DP:GQ"]
  }
  dt
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

if (scenario == "pops") {
  n_pop <- 4
  n_per_pop <- 100
  dir.create(file.path(out_dir, "pops"), recursive = TRUE)
  pop_files <- file.path(out_dir, "pops", sprintf("pop%d.vcf", seq_len(n_pop)))
  for (k in seq_len(n_pop)) {
    writeLines(header(sprintf("pop%d_%03d", k, seq_len(n_per_pop))), pop_files[k])
  }
  for (chr in chroms) {
    sites <- make_sites(chr, round(1.05e6 / 22))
    p <- runif(nrow(sites), 0.01, 0.5)
    for (k in seq_len(n_pop)) {
      in_file <- runif(nrow(sites)) < 0.97
      file_sites <- copy(sites[in_file])
      p_k <- pmin(pmax(p[in_file] + rnorm(sum(in_file), 0, 0.05), 0), 1)
      dosage <- make_dosages(p_k, n_per_pop)
      if (k > 1) {
        # swapped REF/ALT (genotypes recoded), and a different ALT
        swap <- runif(nrow(file_sites)) < 0.01
        file_sites[swap, c("REF", "ALT") := .(ALT, REF)]
        dosage[swap, ] <- 2L - dosage[swap, ]
        other_alt <- !swap & runif(nrow(file_sites)) < 0.001
        bases <- c("A", "C", "G", "T")
        file_sites[other_alt, ALT := vapply(seq_len(.N), function(i) {
          setdiff(bases, c(REF[i], ALT[i]))[1]
        }, character(1))]
      }
      append_rows(pop_files[k], file_sites, gt_strings(dosage))
    }
  }
} else if (scenario %in% c("sgdp", "many")) {
  info <- scenario == "sgdp"
  n_ind <- if (info) 50 else 200
  n_sites <- if (info) 2.5e6 else 1e6
  dir.create(file.path(out_dir, scenario), recursive = TRUE)
  ind_files <- file.path(
    out_dir, scenario, sprintf("S%03d.vcf", seq_len(n_ind))
  )
  for (i in seq_len(n_ind)) {
    writeLines(header(sprintf("S%03d", i), info = info), ind_files[i])
  }
  for (chr in chroms) {
    sites <- make_sites(chr, round(n_sites / 22), info = info)
    dosage <- make_dosages(runif(nrow(sites), 0.01, 0.5), n_ind)
    for (i in seq_len(n_ind)) {
      has_alt <- dosage[, i] > 0
      gs <- c("0/1", "1/1")[dosage[has_alt, i]]
      if (info) {
        n_alt <- sum(has_alt)
        gs <- sprintf("%s:%d:%d", gs, rpois(n_alt, 30), sample.int(99, n_alt, TRUE))
      }
      append_rows(ind_files[i], sites[has_alt], matrix(gs, ncol = 1))
    }
  }
} else {
  stop("unknown scenario ", scenario)
}
cat("done", scenario, "\n")
