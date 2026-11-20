# Time tidypopgen's KING-robust on HGDP650, one task per process.
# Usage: Rscript time_tpg.R <lib> <task> <n_cores> <label> <out.tsv>
# task: full (pairwise_king() matrix), full_kernel (only the C++ calls of
# snp_king(): increment_king_numerator() over the blocks of loci), grouped
# (pairwise_king() within populations, list of matrices), pack (packing
# alone: king_counts_within_groups() with groups of one, so no pairs).
# Each task is run twice in the process: "first" (after gt_load(), as a user
# would) and "second" (warm). If PERF_CTL / PERF_ACK are set, perf counts
# only the first call (perf stat --control fifo:ctl,ack --delay=-1).
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
task <- args[2]
n_cores <- as.integer(args[3])
label <- args[4]
out <- args[5]
suppressPackageStartupMessages(library(tidypopgen, lib.loc = lib))
d <- "/media/andrea/annex/git/tidypopgen/vignettes/articles/data/hgdp"
t_load <- system.time({
  gt <- gt_load(
    "/media/andrea/annex/tmp/claude/tidypopgen/bench_inb/data/hgdp650.gt"
  )
})[["elapsed"]]
pops <- read.delim(file.path(d, "hgdp650_id_pop.txt"))
gt$population <- pops$population[match(gt$id, pops$Id)]
gt_grp <- dplyr::group_by(gt, population)

X <- attr(gt$genotypes, "fbm") # nolint
ind_row <- vctrs::vec_data(gt$genotypes)
ind_col <- show_loci(gt$genotypes)$big_index
n <- length(ind_row)

f <- switch(task,
  full = function() pairwise_king(gt, as_matrix = TRUE, n_cores = n_cores),
  grouped = function() {
    pairwise_king(gt_grp, as_matrix = TRUE, n_cores = n_cores)
  },
  full_kernel = function() {
    # as snp_king(), without the division and the R matrices
    k <- bigstatsr::FBM(n, n, init = 0)
    n_aa <- bigstatsr::FBM(n, n, init = 0)
    intervals <- tidypopgen:::CutBySize(length(ind_col),
      bigstatsr::block_size(n))
    for (j in bigstatsr::rows_along(intervals)) {
      tidypopgen:::increment_king_numerator(k, n_aa, X, ind_row,
        ind_col[seq(intervals[j, 1], intervals[j, 2])], ncores = n_cores)
    }
    k
  },
  pack = function() {
    # as .pairwise_king_by_group() with one individual per group: all the
    # packing, no pairs
    intervals <- tidypopgen:::CutBySize(length(ind_col),
      bigstatsr::block_size(n))
    counts <- 0
    for (j in bigstatsr::rows_along(intervals)) {
      counts <- counts + tidypopgen:::king_counts_within_groups(X, ind_row,
        ind_col[seq(intervals[j, 1], intervals[j, 2])], seq_len(n),
        n_cores)
    }
    counts
  },
  stop("unknown task")
)

ctl <- Sys.getenv("PERF_CTL")
ack <- Sys.getenv("PERF_ACK")
perf_command <- function(cmd) {
  if (!nzchar(ctl)) return(invisible())
  con <- file(ctl, "w", raw = TRUE)
  writeLines(cmd, con)
  close(con)
  con <- file(ack, "r", raw = TRUE)
  readLines(con, n = 1)
  close(con)
}

perf_command("enable")
t1 <- system.time(res <- f())
perf_command("disable")
t2 <- system.time(res <- f())
line <- sprintf("%s\t%s\t%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n",
  format(Sys.time(), "%H:%M:%S"), label, task, n_cores, t_load,
  t1[["elapsed"]], t2[["elapsed"]], t1[["user.self"]] + t1[["user.child"]],
  sub(".*average: ", "", system("uptime", intern = TRUE)))
cat(line)
cat(line, file = out, append = TRUE)
