# Time only the .bed -> .bk step (fbm_read_bed()) of a gen_tibble() reader
# usage: Rscript reader_only.R <lib> <bed file> <out dir> [n_cores]
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(tidypopgen, lib.loc = args[1]))
ns <- asNamespace("tidypopgen")
bed <- args[2]
out <- file.path(args[3], "reader")
n_cores <- if (length(args) > 3) as.integer(args[4]) else 1L
fam <- sub("\\.bed$", ".fam", bed)
n_indiv <- length(readLines(fam))
n_snp <- (file.size(bed) - 3) / ceiling(n_indiv / 4)
stopifnot(n_snp == round(n_snp))
# warm up the bigstatsr reference classes (their first use costs ~0.4 s)
tmp <- bigstatsr::FBM.code256(2, 2, code = bigsnpr::CODE_012,
                              backingfile = tempfile(tmpdir = args[3]))
unlink(paste0(out, c(".bk", ".rds")))
has_n_cores <- "n_cores" %in% names(formals(ns$fbm_read_bed))
t0 <- proc.time()
if (has_n_cores) {
  ns$fbm_read_bed(bed, n_indiv = n_indiv, n_snp = n_snp, backingfile = out,
                  n_cores = n_cores)
} else {
  ns$fbm_read_bed(bed, n_indiv = n_indiv, n_snp = n_snp, backingfile = out)
}
t1 <- proc.time() - t0
cat(sprintf("fbm_read_bed n_indiv %d n_snp %d cores %d: elapsed %.2f user %.2f sys %.2f\n",
            n_indiv, n_snp, if (has_n_cores) n_cores else 1L,
            t1[["elapsed"]], t1[["user.self"]], t1[["sys.self"]]))
if (nzchar(Sys.getenv("BK_MD5"))) {
  cat("md5", unname(tools::md5sum(paste0(out, ".bk"))), "\n")
}
unlink(paste0(out, c(".bk", ".rds")))
unlink(tmp$backingfile)
