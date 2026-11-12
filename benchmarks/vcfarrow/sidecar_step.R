# One measurement of the cost of the locus sidecar (names_as_int = TRUE), in
# this (fresh) R process, for one installed version of tidypopgen. Run by
# sidecar_bench.py, which samples the memory of the process.
#   Rscript sidecar_step.R <lib> <work> <vcf> <threads> <results.jsonl> <action>
# Actions:
#   read_int     gen_tibble(names_as_int = TRUE) into <work>/gt (then
#                gt_save() to <work>/gt.gt, timed apart), with the sizes of
#                the .bk, .rds and .loci.tsv
#   read_chr     the same with names_as_int = FALSE (into <work>/gt_chr, which
#                is deleted), for reference
#   lookup       gt_lookup_loci() of 100 random rows: first call (which builds
#                the sparse index), then warm calls; then warm lookups of
#                1e3..1e6 random rows, and sidecar_index() alone
#   lookup_cold  only the first lookup of 100 random rows (sidecar_bench.py
#                evicts the sidecar from the page cache first)
#   restore      gt_restore_loci() on all loci
#   restore_unchecked  gt_restore_loci(check_positions = FALSE)
#   restore_parts  the parts of gt_restore_loci(): the fetch of all the names
#                and their assignment to the loci table; then the fetch again,
#                with and without the check of positions
#   update       gt_update_backingfile() of 10% of the loci (random), and,
#                with the sidecar, copy_loci_sidecar() alone
# Every action but read_* starts from gt_load(<work>/gt.gt).
# Before each timed step, "<step>\t<epoch>" is written to <work>/current_step
# (as in run_step.R).

args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
work <- args[2]
vcf <- args[3]
threads <- as.integer(args[4])
results <- args[5]
action <- args[6]

suppressPackageStartupMessages(library(tidypopgen, lib.loc = lib))
suppressPackageStartupMessages(library(dplyr))
version <- as.character(utils::packageVersion("tidypopgen", lib.loc = lib))
ns <- asNamespace("tidypopgen")
has_sidecar <- exists("gt_lookup_loci", envir = ns)
loci_header <- "int_name\tname\tchromosome\tposition"

record <- function(step, t0, t1, tm, status, summary, extra = list()) {
  line <- jsonlite::toJSON(list(
    step = step, phase = "step", t_start = t0, t_end = t1,
    # (wall time from Sys.time(), which resolves microseconds; system.time()
    # only resolves milliseconds)
    elapsed = t1 - t0,
    user = unname(tm[["user.self"]] + tm[["user.child"]]),
    sys = unname(tm[["sys.self"]] + tm[["sys.child"]]),
    status = status, summary = summary, version = version,
    extra = extra
  ), auto_unbox = TRUE, digits = NA, na = "null")
  cat(line, "\n", file = results, append = TRUE, sep = "")
}

# time fun(), which returns list(value, summary, extra); returns the value
# (NULL on error)
timed <- function(step, fun, gc_first = TRUE) {
  writeLines(sprintf("%s\t%.3f", step, as.numeric(Sys.time())),
    file.path(work, "current_step"))
  if (gc_first) gc(full = TRUE)
  t0 <- as.numeric(Sys.time())
  res <- NULL
  status <- "ok"
  err <- NA_character_
  tm <- system.time(
    res <- tryCatch(fun(), error = function(e) {
      status <<- "error"
      err <<- conditionMessage(e)
      NULL
    }),
    # (gc() was run above, if wanted, outside of [t0, t1])
    gcFirst = FALSE
  )
  t1 <- as.numeric(Sys.time())
  summary <- if (status == "ok") res$summary else err
  extra <- if (status == "ok" && !is.null(res$extra)) res$extra else list()
  record(step, t0, t1, tm, status, summary, extra)
  if (status != "ok") {
    return(NULL)
  }
  res$value
}

fsize <- function(path) if (file.exists(path)) file.size(path) else NA_real_
gt_file <- file.path(work, "gt.gt")

# a short check of looked-up names: the first ones and their number
names_summary <- function(nm) {
  sprintf("%d names (%d NA), first: %s", length(nm), sum(is.na(nm)),
    paste(utils::head(nm, 3), collapse = ", "))
}

if (action %in% c("read_int", "read_chr")) {
  int <- action == "read_int"
  bk <- file.path(work, if (int) "gt" else "gt_chr")
  unlink(paste0(bk, c(".bk", ".rds", ".gt", ".loci.tsv")))
  x <- timed(action, function() {
    x <- gen_tibble(vcf,
      backingfile = bk, parser = "cpp", n_cores = threads,
      names_as_int = int, quiet = TRUE
    )
    list(
      value = x,
      summary = sprintf("%d indiv x %d loci", nrow(x), count_loci(x)),
      extra = list(
        n_loci = count_loci(x),
        bk_bytes = fsize(paste0(bk, ".bk")),
        rds_bytes = fsize(paste0(bk, ".rds")),
        sidecar_bytes = fsize(paste0(bk, ".loci.tsv"))
      )
    )
  })
  if (is.null(x)) quit(status = 1)
  if (int) {
    # kept for the other actions (and not part of the read time)
    timed("gt_save", function() {
      gt_save(x, file_name = gt_file, quiet = TRUE)
      list(value = NULL, summary = "saved",
        extra = list(gt_bytes = fsize(gt_file)))
    })
  } else {
    rm(x)
    gc()
    unlink(paste0(bk, c(".bk", ".rds", ".gt")))
  }
  quit(status = 0)
}

x <- timed("gt_load", function() {
  list(value = gt_load(gt_file), summary = "loaded")
})
if (is.null(x)) quit(status = 1)
n <- count_loci(x)

# warm calls of fun(rows), with new random rows each time: their wall times
# (from Sys.time(), to the microsecond)
warm_times <- function(fun, size, reps, seed) {
  set.seed(seed)
  vapply(seq_len(reps), function(i) {
    rows <- sample.int(n, size)
    t0 <- as.numeric(Sys.time())
    fun(rows)
    as.numeric(Sys.time()) - t0
  }, numeric(1))
}

if (action %in% c("lookup", "lookup_cold")) {
  if (!has_sidecar) stop("this version has no gt_lookup_loci()")
  set.seed(1)
  rows1 <- sample.int(n, 100)
  timed("lookup100_first", function() {
    res <- gt_lookup_loci(x, rows1)
    list(value = res, summary = names_summary(res$name),
      extra = list(index_entries = length(ns$sidecar_cache[[
        normalizePath(ns$loci_sidecar_path(x))]]$index$keys)))
  })
  if (action == "lookup_cold") quit(status = 0)
  # the same 100 rows again (index cached, lines in the page cache)
  timed("lookup100_same", function() {
    res <- gt_lookup_loci(x, rows1)
    list(value = res, summary = names_summary(res$name))
  })
  # new random rows, index cached: 20 calls, each timed
  for (size in c(100L, 1000L, 10000L, 100000L, 1000000L)) {
    if (size > n) next
    reps <- if (size <= 10000L) 20L else 5L
    timed(sprintf("lookup%d_warm", size), function() {
      t <- warm_times(function(rows) gt_lookup_loci(x, rows), size, reps,
        seed = 10 + size)
      list(value = NULL,
        summary = sprintf("%d calls: median %.4f s", reps, median(t)),
        extra = list(times = t, median = median(t), min = min(t),
          max = max(t)))
    }, gc_first = FALSE)
    # the C++ lookup alone (sorted integer keys, positions checked), to
    # separate it from the R subsetting of the loci table
    keys_all <- show_loci(x)$name
    pos_all <- show_loci(x)$position
    timed(sprintf("fetch%d_warm", size), function() {
      t <- warm_times(function(rows) {
        ns$loci_sidecar_names(x$genotypes, keys_all[rows], pos_all[rows])
      }, size, reps, seed = 10 + size)
      list(value = NULL,
        summary = sprintf("%d calls: median %.4f s", reps, median(t)),
        extra = list(times = t, median = median(t), min = min(t),
          max = max(t)))
    }, gc_first = FALSE)
  }
  # building the sparse index alone (no cache)
  path <- ns$loci_sidecar_path(x)
  timed("sidecar_index", function() {
    t <- vapply(1:3, function(i) {
      t0 <- as.numeric(Sys.time())
      ns$sidecar_index(path, loci_header, 1024L)
      as.numeric(Sys.time()) - t0
    }, numeric(1))
    list(value = NULL, summary = sprintf("median %.3f s", median(t)),
      extra = list(times = t, median = median(t),
        sidecar_bytes = fsize(path)))
  })
  quit(status = 0)
}

if (action == "restore") {
  if (!has_sidecar) stop("this version has no gt_restore_loci()")
  timed("restore", function() {
    y <- gt_restore_loci(x)
    list(value = y, summary = names_summary(show_loci(y)$name))
  })
  quit(status = 0)
}

if (action == "restore_unchecked") {
  if (!has_sidecar) stop("this version has no gt_restore_loci()")
  timed("restore_unchecked", function() {
    y <- gt_restore_loci(x, check_positions = FALSE)
    list(value = y, summary = names_summary(show_loci(y)$name))
  })
  quit(status = 0)
}

if (action == "restore_parts") {
  # where the time of gt_restore_loci() goes: the C++ fetch of all the names
  # (with the check of the positions), then the assignment to the loci table;
  # then the fetch again (once R has made room for the strings), with and
  # without the check
  if (!has_sidecar) stop("this version has no gt_restore_loci()")
  keys <- show_loci(x)$name
  pos <- show_loci(x)$position
  nm <- timed("restore_fetch_checked_1", function() {
    nm <- ns$loci_sidecar_names(x$genotypes, keys, pos)
    list(value = nm, summary = names_summary(nm))
  })
  timed("restore_assign", function() {
    show_loci(x)$name <- nm
    list(value = NULL, summary = names_summary(show_loci(x)$name))
  })
  rm(nm)
  nm <- timed("restore_fetch_checked_2", function() {
    nm <- ns$loci_sidecar_names(x$genotypes, keys, pos)
    list(value = nm, summary = names_summary(nm))
  })
  rm(nm)
  nm <- timed("restore_fetch_unchecked_3", function() {
    nm <- ns$loci_sidecar_names(x$genotypes, keys, NULL)
    list(value = nm, summary = names_summary(nm))
  })
  quit(status = 0)
}

if (action == "update") {
  out <- file.path(work, "out")
  dir.create(out, showWarnings = FALSE)
  unlink(list.files(out, full.names = TRUE))
  set.seed(3)
  keep <- sort(sample.int(n, round(n / 10)))
  x10 <- timed("select10", function() {
    y <- select_loci(x, .sel_arg = all_of(keep))
    list(value = y, summary = sprintf("%d loci", count_loci(y)))
  })
  if (is.null(x10)) quit(status = 1)
  new_bk <- file.path(out, "sub")
  timed("update10", function() {
    y <- gt_update_backingfile(x10, backingfile = new_bk, quiet = TRUE)
    list(value = NULL, summary = sprintf("%d loci", count_loci(y)),
      extra = list(bk_bytes = fsize(paste0(new_bk, ".bk")),
        sidecar_bytes = fsize(paste0(new_bk, ".loci.tsv"))))
  })
  if (has_sidecar) {
    # the copy of the sidecar alone (part of update10 above)
    timed("sidecar_copy10", function() {
      bk2 <- file.path(out, "sub2.bk")
      written <- ns$copy_loci_sidecar(x10, bk2, show_loci(x10)$name)
      path2 <- ns$sidecar_path_from_bk(bk2, "loci")
      list(value = NULL, summary = sprintf("written: %s", written),
        extra = list(sidecar_bytes = fsize(path2)))
    })
    # check the new sidecar: the names of the first and last loci kept
    timed("check_update10", function() {
      y <- gt_load(file.path(out, "sub.gt"))
      a <- gt_lookup_loci(y, c(1, count_loci(y)))$name
      b <- gt_lookup_loci(x, keep[c(1, length(keep))])$name
      list(value = NULL, summary = sprintf("%s; identical: %s",
        paste(a, collapse = ", "), identical(a, b)))
    })
  }
  unlink(list.files(out, full.names = TRUE))
  quit(status = 0)
}

stop("unknown action ", action)
