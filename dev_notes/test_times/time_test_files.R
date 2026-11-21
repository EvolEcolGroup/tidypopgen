# Time each test file of tidypopgen against an installed (-O2) build, as
# CRAN runs them: NOT_CRAN unset (skip_on_cran() applies), non-interactive.
#
# Environment variables:
#   TT_LIB  library holding the -O2 build of tidypopgen
#   TT_SRC  source tree whose tests/testthat is run
#   TT_OUT  folder for per_file.csv and per_test.csv

lib <- Sys.getenv("TT_LIB")
src <- Sys.getenv("TT_SRC")
out <- Sys.getenv("TT_OUT")

# use the -O2 build first
.libPaths(c(lib, .libPaths()))
pkg_path <- find.package("tidypopgen")
stopifnot(normalizePath(dirname(pkg_path)) == normalizePath(lib))
stopifnot(Sys.getenv("NOT_CRAN") == "")
cat("tidypopgen", format(utils::packageVersion("tidypopgen")), "from", pkg_path, "\n")

test_dir <- file.path(src, "tests", "testthat")
files <- sort(list.files(test_dir, pattern = "^test.*\\.R$"))

per_file <- list()
per_test <- list()
for (f in files) {
  gc()
  # time the whole file, including code outside test_that() blocks
  timing <- system.time(
    results <- tryCatch(
      testthat::test_file(
        file.path(test_dir, f),
        reporter = "silent",
        package = "tidypopgen",
        load_package = "installed"
      ),
      error = function(e) {
        message("error in ", f, ": ", conditionMessage(e))
        NULL
      }
    )
  )
  tests <- if (is.null(results)) NULL else as.data.frame(results)
  per_file[[f]] <- data.frame(
    file = f,
    elapsed_s = unname(timing[["elapsed"]]),
    cpu_s = unname(timing[["user.self"]] + timing[["sys.self"]] +
      timing[["user.child"]] + timing[["sys.child"]]),
    n_tests = if (is.null(tests)) NA else nrow(tests),
    n_skipped = if (is.null(tests)) NA else sum(tests$skipped),
    n_failed = if (is.null(tests)) NA else sum(tests$failed > 0 | tests$error)
  )
  if (!is.null(tests)) {
    per_test[[f]] <- tests[, c("file", "test", "real", "user", "system",
                               "skipped", "failed", "error")]
  }
  cat(sprintf("%-50s %7.1f s\n", f, timing[["elapsed"]]))
  # write as we go, so a crash keeps what was timed
  utils::write.csv(do.call(rbind, per_file), file.path(out, "per_file.csv"),
                   row.names = FALSE)
}
utils::write.csv(do.call(rbind, per_test), file.path(out, "per_test.csv"),
                 row.names = FALSE)
cat("total elapsed", sum(sapply(per_file, `[[`, "elapsed_s")), "s\n")
