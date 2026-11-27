# Rscript time_proto.R <bgzipped file>: time ways of decompressing it
args <- commandArgs(TRUE)
proto_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
Sys.setenv(R_MAKEVARS_USER = file.path(proto_dir, "Makevars_proto"))
Rcpp::sourceCpp(file.path(proto_dir, "bgzf_proto.cpp"), cacheDir = file.path(proto_dir, "cache"))
f <- args[1]
tm <- function(label, expr) {
  expr <- substitute(expr)
  for (rep in 1:3) {
    t <- system.time(r <- eval(expr))
    cat(sprintf("%-26s rep %d  %6.2f s elapsed  %6.2f s user  (%0.f bytes, %0.f lines)\n", label, rep, t[["elapsed"]], t[["user.self"]], r[1], r[2]))
  }
}
tm("zlib stream (gzread)", zlib_stream(f))
for (n in c(1, 2, 4)) tm(sprintf("bgzf zlib, %d threads", n), bgzf_parallel(f, n, FALSE))
for (n in c(1, 2, 4)) tm(sprintf("bgzf libdeflate, %d threads", n), bgzf_parallel(f, n, TRUE))
for (lib in c(FALSE, TRUE)) for (n in c(1, 4)) cat(sprintf("crc32 %-10s %d threads: %.2f s\n", if (lib) "libdeflate" else "zlib", n, crc_cost(f, n, lib)))
