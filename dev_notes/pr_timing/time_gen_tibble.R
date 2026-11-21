# Time gen_tibble() ("cpp" parser) on one VCF, in a fresh R process.
# Usage: Rscript time_gen_tibble.R <lib> <vcf> <n_cores> <label> <out.tsv>
# Appends: date, label, commit, file, n_cores, wall s, user s, load (1 min),
# md5 of the .bk (to check that the versions give the same genotypes)
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
vcf <- args[2]
n_cores <- as.integer(args[3])
label <- args[4]
out <- args[5]
suppressPackageStartupMessages(library(tidypopgen, lib.loc = lib))
commit_file <- file.path(lib, "COMMIT")
commit <- if (file.exists(commit_file)) readLines(commit_file)[1] else NA
load1 <- strsplit(readLines("/proc/loadavg"), " ")[[1]][1]
bk <- tempfile("gt_")
timing <- system.time(
  gt <- gen_tibble(vcf,
    parser = "cpp", backingfile = bk, n_cores = n_cores, quiet = TRUE
  )
)
md5 <- unname(substr(tools::md5sum(paste0(bk, ".bk")), 1, 8))
row <- data.frame(
  date = format(Sys.time(), "%FT%T"), version = label, commit = commit,
  file = basename(vcf), n_cores = n_cores,
  wall_s = round(timing[["elapsed"]], 3),
  user_s = round(timing[["user.self"]] + timing[["user.child"]], 3),
  load1 = load1, bk_md5 = md5
)
utils::write.table(row, out,
  sep = "\t", quote = FALSE, row.names = FALSE,
  col.names = !file.exists(out), append = file.exists(out)
)
unlink(paste0(bk, c(".bk", ".rds", ".gt", ".loci.tsv")))
cat(label, basename(vcf), n_cores, timing[["elapsed"]], md5, "\n")
