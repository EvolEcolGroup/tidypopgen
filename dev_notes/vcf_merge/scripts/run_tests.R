# Rscript run_tests.R <lib> <test files...>: run test files of the worktree
# against an installed tidypopgen
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
wt <- file.path(Sys.getenv("TPG_WORK", "/media/andrea/Elements1/tmp/claude/tidypopgen"), "merge_vcf/wt/tests/testthat")
res <- lapply(args[-1], function(f) {
  testthat::test_file(file.path(wt, f), package = "tidypopgen",
    load_package = "installed", reporter = "summary")
})
df <- do.call(rbind, lapply(res, as.data.frame))
cat("failed:", sum(df$failed), " errors:", sum(df$error), " skipped:",
    sum(df$skipped), " passed:", sum(df$nb) - sum(df$failed), "\n")
