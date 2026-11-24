#!/bin/bash
# usage: run_tests.sh <source dir> <library> <output prefix> [test filter]
# installs the package from <source dir> into <library> and runs its tests
set -euo pipefail
src=$1; lib=$2; out=$3; filter=${4:-}
mkdir -p "$lib"
R CMD INSTALL --no-docs --no-test-load --library="$lib" "$src" > "$out.install.log" 2>&1
cd "$src/tests/testthat"
R_LIBS="$lib" Rscript -e '
args <- commandArgs(trailingOnly = TRUE)
library(testthat)
res <- as.data.frame(test_dir(".", package = "tidypopgen", load_package = "installed",
  filter = if (nzchar(args[2])) args[2] else NULL,
  reporter = "summary", stop_on_failure = FALSE))
res$status <- ifelse(res$error | res$failed > 0, "FAIL", ifelse(res$skipped, "skip", "ok"))
write.table(res[, c("file", "test", "nb", "failed", "skipped", "error", "status")],
  paste0(args[1], ".results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(res$status))
' "$out" "$filter" > "$out.tests.log" 2>&1
