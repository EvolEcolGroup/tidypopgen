# Rscript bench_one.R <lib> <geno file> <n_cores> <mode: reader|gen_tibble> <out .bk prefix> [max_buffer_bytes]
args <- commandArgs(TRUE)
lib <- args[1]; file <- args[2]; n_cores <- as.integer(args[3]); mode <- args[4]; bk <- args[5]
max_buf <- if (length(args) >= 6) as.numeric(args[6]) else 2^31
suppressMessages(library(tidypopgen, lib.loc = lib))
n <- 4866; m <- 1233013
is_new <- exists("fbm_read_packedancestry", envir = asNamespace("tidypopgen"))
if (mode == "reader") {
  t <- system.time({
    if (is_new) {
      tidypopgen:::fbm_read_packedancestry(file, n_indiv = n, n_snp = m,
        backingfile = bk, n_cores = n_cores, max_buffer_bytes = max_buf)
    } else {
      fbm <- bigstatsr::FBM.code256(nrow = n, ncol = m, code = bigsnpr::CODE_012, backingfile = bk)
      ok <- tidypopgen:::read_packedancestry(file, fbm, tab = tidypopgen:::get_packedancestry_code())
    }
  })
} else {
  t <- system.time(gt <- gen_tibble(file, backingfile = bk, n_cores = n_cores, quiet = TRUE))
}
cat(sprintf("RESULT lib=%s file=%s cores=%d mode=%s buf=%g elapsed=%.2f user=%.2f sys=%.2f\n",
  basename(lib), basename(file), n_cores, mode, max_buf, t[["elapsed"]], t[["user.self"]], t[["sys.self"]]))
