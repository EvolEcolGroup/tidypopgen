# Run steps of steps.R in this (fresh) R process, for one installed version.
#   Rscript run_step.R <lib> <work> <vcf> <threads> <results.jsonl> <step>...
# read_vcf (followed by set_groups) creates <work>/gt.gt; every other step
# starts from gt_load(<work>/gt.gt), so each step is timed (and its memory
# measured by bench.py) on its own, from the unfiltered data.
# Before each step, "<step>\t<epoch>" is written to <work>/current_step, so
# that bench.py can enforce the time limit per step.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
work <- args[2]
vcf <- args[3]
threads <- as.integer(args[4])
results <- args[5]
step_ids <- args[-(1:5)]

suppressPackageStartupMessages(library(tidypopgen, lib.loc = lib))
here <- dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)
)))
source(file.path(here, "steps.R"))

out_dir <- file.path(work, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ctx <- list(
  threads = threads, out_dir = out_dir, work = work, vcf = vcf,
  strata = sub("\\.vcf(\\.gz)?$", ".strata.tsv", vcf)
)
gt_file <- file.path(work, "gt.gt")

record <- function(step, phase, t0, t1, tm, status, summary) {
  line <- jsonlite::toJSON(list(
    step = step, phase = phase, t_start = t0, t_end = t1,
    elapsed = unname(tm[["elapsed"]]),
    user = unname(tm[["user.self"]] + tm[["user.child"]]),
    sys = unname(tm[["sys.self"]] + tm[["sys.child"]]),
    status = status, summary = summary,
    version = as.character(utils::packageVersion("tidypopgen", lib.loc = lib))
  ), auto_unbox = TRUE, digits = NA, na = "null")
  cat(line, "\n", file = results, append = TRUE, sep = "")
}

# runs fun(), timing it; returns list(value, ok)
timed <- function(step, phase, fun) {
  writeLines(sprintf("%s\t%.3f", step, as.numeric(Sys.time())),
    file.path(work, "current_step"))
  gc(full = TRUE)
  t0 <- as.numeric(Sys.time())
  value <- NULL
  status <- "ok"
  summary <- NA_character_
  tm <- system.time(
    value <- tryCatch(fun(), error = function(e) {
      status <<- "error"
      summary <<- conditionMessage(e)
      NULL
    })
  )
  t1 <- as.numeric(Sys.time())
  if (status == "ok") {
    summary <- if (is.character(value)) value else if (inherits(value, "gen_tbl")) {
      dims(value)
    } else {
      class(value)[1]
    }
  }
  record(step, phase, t0, t1, tm, status, summary)
  list(value = value, ok = status == "ok")
}

if ("read_vcf" %in% step_ids) {
  unlink(list.files(work, pattern = "^gt\\.", full.names = TRUE))
  res <- timed("read_vcf", "step", function() steps$read_vcf$fun(NULL, ctx))
  if (!res$ok) quit(status = 1)
  res <- timed("set_groups", "step", function() steps$set_groups$fun(res$value, ctx))
  if (!res$ok) quit(status = 1)
  saved <- timed("gt_save", "io", function() {
    gt_save(res$value, file_name = gt_file, quiet = TRUE)
    "saved"
  })
  if (!saved$ok) quit(status = 1)
  step_ids <- setdiff(step_ids, c("read_vcf", "set_groups"))
}
if (length(step_ids) == 0) quit(status = 0)

loaded <- timed("gt_load", "io", function() gt_load(gt_file))
if (!loaded$ok) quit(status = 1)
gt <- loaded$value

for (s in step_ids) {
  st <- steps[[s]]
  if (is.null(st)) stop("unknown step ", s)
  if (!is.null(st$needs)) {
    gt_args <- names(formals(utils::getS3method("gen_tibble", "character",
      envir = asNamespace("tidypopgen"))))
    if (!all(st$needs %in% gt_args)) {
      record(s, "step", NA, NA, c(elapsed = NA, user.self = NA,
        user.child = NA, sys.self = NA, sys.child = NA),
        "n/a", "not available in this version")
      next
    }
  }
  timed(s, "step", function() st$fun(gt, ctx))
}
