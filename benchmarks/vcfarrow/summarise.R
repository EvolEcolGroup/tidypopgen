# Summarise the results of bench.py: old vs new tidypopgen on the same data
# and machine, and an extrapolation to the full dataset (18,579,683 variants)
# next to the numbers published in the VCFArrow comparison (4 threads /
# 4 workers).
#   Rscript summarise.R [results/*.tsv] > results/summary.md

suppressPackageStartupMessages(library(dplyr))
here <- dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)
)))
files <- commandArgs(trailingOnly = TRUE)
if (length(files) == 0) {
  files <- list.files(file.path(here, "results"), "\\.tsv$", full.names = TRUE)
}
res <- bind_rows(lapply(files, function(f) {
  utils::read.delim(f, na.strings = "NA", colClasses = c(summary = "character"))
}))
full_n <- 18579683

# published: tidypopgen (4 threads) and VCFArrow (4 workers); time in s,
# memory median/peak PSS in GB (one session for all functions, so the
# memory includes everything loaded before)
pub <- tibble::tribble(
  ~step, ~vcfarrow, ~tpg_s, ~tpg_mem, ~va_s, ~va_mem,
  "read_vcf", "read VCF", 30 * 60, "23.38/32.29", 6.9 * 60, "8.94/14.10",
  "set_groups", "set groups", 0, "30.96/30.96", 0, "13.06/13.06",
  "extract_samples", "extract samples", 0, "30.96/30.96", 55.6, "16.65/19.05",
  "extract_groups", "extract groups", 0, "30.96/30.96", 60, "19.20/21.61",
  "subsample_loci", "subsample loci (1000)", 6.1 * 60, "30.96/30.96", 0.5, "21.06/21.11",
  "filter_maf", "filter MAF (0.05)", 24.1, "30.96/38.25", 33.7, "18.47/23.73",
  "filter_locus_missing", "filter locus missing (0.2)", 21.0, "30.96/38.25", 19.0, "18.27/23.94",
  "filter_sample_missing", "filter sample missing (0.8)", 4.7, "30.96/30.96", 3.6 * 60, "20.44/24.04",
  "filter_invariant", "filter invariant", 20.1, "30.96/38.25", 41.8, "19.63/24.21",
  "merge", "merge two VCFs", 1.4 * 60, "32.65/38.25", 10.3 * 60, "31.79/32.53",
  "sample_stats", "sample statistics", 4.7, "30.96/33.82", 60, "23.35/32.27",
  "theta_pi", "theta / pi", 0.6, "30.96/30.96", 1.1 * 60, "24.21/25.10",
  "missing_plot", "missing data plot + table", 4.7, "30.96/30.96", 46.0, "13.64/13.64",
  "write_vcf", "write VCF", Inf, "30.96/30.96", 8.5 * 60, "15.87/16.64",
  "to_genlight", "-> genlight", 3.3 * 60, "38.25/39.01", 6.9 * 60, "36.18/42.33",
  "to_ped", "-> PLINK ped", Inf, "32.53/38.64", 2.9 * 60, "29.10/41.02",
  "to_bed", "-> PLINK bed", 4.1, "30.84/31.10", 1.7 * 60, "27.28/41.02"
)
extra_steps <- c("subsample_loci_bare", "theta_pi_grouped",
  "filter_chain_update", "read_vcf_int")
step_order <- c(pub$step, extra_steps)

fmt_t <- function(s) {
  ifelse(is.na(s), "—", ifelse(is.infinite(s), "> 60 min",
    ifelse(s < 60, sprintf("%.1f s", s), ifelse(s < 3600,
      sprintf("%.1f min", s / 60), sprintf("%.2f h", s / 3600)))))
}
fmt_gb <- function(x) ifelse(is.na(x), "—", sprintf("%.2f", x))
md_table <- function(df) {
  cat("| ", paste(names(df), collapse = " | "), " |\n", sep = "")
  cat("|", paste(rep("---", ncol(df)), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(df))) {
    cat("| ", paste(unlist(df[i, ]), collapse = " | "), " |\n", sep = "")
  }
  cat("\n")
}

# versions in the order they were run (e.g. 0.4.4, 1.0.0, then a branch),
# before group_by() sorts them alphabetically: later versions are compared
# with the first
labels <- unique(res$label)
res <- res %>%
  # a step that was run again replaces the earlier run
  group_by(label, n_variants, threads, step) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  mutate(
    elapsed_s = ifelse(grepl("^timeout", status), Inf, elapsed_s),
    step = factor(step, levels = unique(c(step_order, step)))
  )

cat("# tidypopgen on the VCFArrow benchmark\n\n")
cat("Synthetic DiscoSnpRad VCF (56 samples), one fresh R process per step,",
  "peak PSS of the process tree (data loaded + step).\n\n")

# 1. same data, same machine -------------------------------------------------
for (n in sort(unique(res$n_variants))) {
  for (th in sort(unique(res$threads))) {
    r <- res %>% filter(n_variants == n, threads == th)
    if (nrow(r) == 0) next
    cat(sprintf("## %s variants, %d threads\n\n", format(n, big.mark = ","), th))
    # the versions run at this size (e.g. 0.4.4 is not run at full size)
    labels_n <- labels[labels %in% r$label]
    wide <- lapply(labels_n, function(l) {
      r %>%
        filter(label == l) %>%
        transmute(step,
          !!paste(l, "time") := ifelse(status == "ok", fmt_t(elapsed_s),
            status),
          !!paste(l, "peak GB") := fmt_gb(pss_peak_gb),
          !!paste0(".t_", l) := elapsed_s,
          !!paste0(".m_", l) := pss_peak_gb,
          !!paste0(".s_", l) := summary)
    })
    tab <- Reduce(function(a, b) full_join(a, b, by = "step"), wide) %>%
      arrange(step)
    # each later version against the first one run at this size (e.g. 0.4.4);
    # with more than two versions, the columns are named after the later one
    ref <- labels_n[1]
    diffs <- character(0)
    for (l in labels_n[-1]) {
      sfx <- if (length(labels_n) > 2) paste0(" ", l) else ""
      t_old <- tab[[paste0(".t_", ref)]]
      t_new <- tab[[paste0(".t_", l)]]
      m_old <- tab[[paste0(".m_", ref)]]
      m_new <- tab[[paste0(".m_", l)]]
      tab[[paste0("speed_up", sfx)]] <- ifelse(
        is.finite(t_old) & is.finite(t_new) & t_new > 0,
        sprintf("%.1fx", t_old / t_new), "")
      tab[[paste0("memory", sfx)]] <- ifelse(!is.na(m_old) & !is.na(m_new),
        sprintf("%.2fx", m_new / m_old), "")
      s_old <- tab[[paste0(".s_", ref)]]
      s_new <- tab[[paste0(".s_", l)]]
      na_old <- is.na(s_old) | grepl("not available", s_old)
      na_new <- is.na(s_new) | grepl("not available", s_new)
      same <- ifelse(na_old | na_new, "",
        ifelse(s_old == s_new, "yes", "**no**"))
      tab[[paste0("same_result", sfx)]] <- same
      differ <- which(same %in% "**no**")
      diffs <- c(diffs, sprintf("- %s (%s / %s): %s / %s",
        as.character(tab$step[differ]), ref, l, s_old[differ], s_new[differ]))
    }
    tab$step <- as.character(tab$step)
    md_table(tab %>% select(-starts_with(".")))
    if (length(diffs) > 0) {
      cat("Results that differ between versions:\n\n")
      cat(diffs, sep = "\n")
      cat("\n\n")
    }
  }
}

# 2. extrapolated to the full dataset ------------------------------------------
# measured at full size where it was run; otherwise, from the two largest
# sizes run (of at least 1e6 variants, so that fixed costs do not dominate),
# a power law t = a * n^b for time and a straight line for peak memory;
# steps that hit the time or memory limit at a smaller size are reported so
extrap <- res %>%
  group_by(label, threads, step) %>%
  summarise(
    t_full = {
      ok <- status == "ok" & is.finite(elapsed_s) & elapsed_s > 0 &
        n_variants >= 1e6
      n_ok <- n_variants[ok]
      if (any(n_variants == full_n & status == "ok")) {
        elapsed_s[n_variants == full_n & status == "ok"][1]
      } else if (any(grepl("^timeout", status))) {
        Inf
      } else if (length(unique(n_ok)) >= 2) {
        keep <- n_ok %in% sort(unique(n_ok), decreasing = TRUE)[1:2]
        fit <- stats::lm(log(elapsed_s[ok][keep]) ~ log(n_ok[keep]))
        exp(stats::predict(fit, data.frame(n_ok = full_n)))
      } else {
        NA_real_
      }
    },
    m_full = {
      ok <- status == "ok" & !is.na(pss_peak_gb) & n_variants >= 1e6
      n_ok <- n_variants[ok]
      if (any(n_variants == full_n & status == "ok")) {
        pss_peak_gb[n_variants == full_n & status == "ok"][1]
      } else if (length(unique(n_ok)) >= 2) {
        keep <- n_ok %in% sort(unique(n_ok), decreasing = TRUE)[1:2]
        fit <- stats::lm(pss_peak_gb[ok][keep] ~ n_ok[keep])
        stats::predict(fit, data.frame(n_ok = full_n))
      } else {
        NA_real_
      }
    },
    status_full = {
      st <- status[n_variants == max(n_variants)][1]
      if (st %in% c("ok", "n/a")) "" else st
    },
    .groups = "drop"
  )
if (nrow(extrap) > 0) {
  for (th in sort(unique(extrap$threads))) {
    cat(sprintf("## Estimated at 18.6 M variants (%d threads), with the published numbers\n\n", th))
    cat("Measured on this machine (4 cores, 15 GB) at",
      paste(format(sort(unique(res$n_variants)), big.mark = ","),
        collapse = ", "), "variants; values at 18.6 M are measured where",
      "that size was run, otherwise extrapolated from the two largest sizes",
      "(\"—\": not enough sizes). The published runs used another machine,",
      "50 GB, a 60-min limit, and one R session for all functions.\n\n")
    tab <- pub %>% select(step, vcfarrow, tpg_s, tpg_mem, va_s)
    for (l in labels) {
      e <- extrap %>% filter(label == l, threads == th) %>%
        mutate(t_full = ifelse(status_full != "" & !is.finite(t_full),
          NA, t_full)) %>%
        select(step, t_full, m_full) %>% mutate(step = as.character(step))
      tab <- tab %>% left_join(e, by = "step") %>%
        rename(!!paste(l, "est. time") := t_full, !!paste(l, "est. peak GB") := m_full)
    }
    out <- tab %>% mutate(
      across(c(tpg_s, va_s, ends_with("est. time")), fmt_t),
      across(ends_with("est. peak GB"), fmt_gb)
    ) %>%
      rename(`published tidypopgen (4)` = tpg_s,
        `published tidypopgen mem` = tpg_mem,
        `published VCFArrow (4)` = va_s, function_ = vcfarrow) %>%
      select(-step)
    names(out)[names(out) == "function_"] <- "function"
    md_table(out)
  }
}
