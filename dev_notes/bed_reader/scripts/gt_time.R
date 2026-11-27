# Time gt_load() of a .gt file, then saving it again (save_gt_file())
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(tidypopgen, lib.loc = args[1]))
t0 <- proc.time()
x <- gt_load(args[2])
t_load <- (proc.time() - t0)[["elapsed"]]
out <- tempfile(tmpdir = args[3], fileext = ".gt")
t0 <- proc.time()
asNamespace("tidypopgen")$save_gt_file(x, out)
t_save <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("load %.2f save %.2f size %.3f GB\n", t_load, t_save,
            file.size(out) / 1e9))
unlink(out)
