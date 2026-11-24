# Ploidy audit: call every analysis function of tidypopgen on simulated
# diploid, haploid, tetraploid, hexaploid and mixed (2x + 4x) data, and record
# whether it runs, warns or errors. Writes audit.tsv next to this script.
#
# usage: Rscript audit.R <library with tidypopgen> <output tsv>
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
out_file <- args[2]
suppressPackageStartupMessages({
  library(tidypopgen, lib.loc = lib)
  library(dplyr)
})
options(bigstatsr.check.parallel.blas = FALSE)
set.seed(301)

n_pop <- 3
n_per_pop <- 12
n_ind <- n_pop * n_per_pop
n_loci <- 600
pop <- rep(paste0("pop", seq_len(n_pop)), each = n_per_pop)
# allele frequencies that differ between populations
p_anc <- runif(n_loci, 0.1, 0.9)
p_pop <- sapply(seq_len(n_pop), function(i) {
  pmin(pmax(p_anc + rnorm(n_loci, 0, 0.12), 0.15), 0.85)
})

sim_gt <- function(ploidy, missing = 0.02) {
  ploidy <- rep_len(ploidy, n_ind)
  geno <- t(sapply(seq_len(n_ind), function(i) {
    rbinom(n_loci, ploidy[i], p_pop[, match(pop[i], unique(pop))])
  }))
  geno[sample(length(geno), length(geno) * missing)] <- NA
  loci <- data.frame(
    name = paste0("rs", seq_len(n_loci)),
    chromosome = rep(1:3, each = n_loci / 3),
    position = rep(seq_len(n_loci / 3) * 5000L, 3),
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
  meta <- data.frame(id = paste0("ind", seq_len(n_ind)), population = pop)
  gen_tibble(
    geno,
    indiv_meta = meta, loci = loci,
    ploidy = if (length(unique(ploidy)) == 1) ploidy[1] else ploidy,
    backingfile = tempfile("audit_"), quiet = TRUE
  )
}

ploidies <- list(
  diploid = 2,
  haploid = 1,
  tetraploid = 4,
  hexaploid = 6,
  mixed_2x_4x = rep(c(2, 4), n_ind / 2)
)
datasets <- lapply(ploidies, sim_gt)
# the same without missing genotypes, for the analyses that need complete
# data (so that they can be audited separately from imputation)
datasets_complete <- lapply(ploidies, sim_gt, missing = 0)

# every call takes the ungrouped gen_tibble x, the grouped one g, and the
# ungrouped one without missing data xc
calls <- list(
  # accessors and basic summaries
  show_ploidy = function(x, g, xc) show_ploidy(x),
  indiv_ploidy = function(x, g, xc) indiv_ploidy(x),
  show_genotypes = function(x, g, xc) show_genotypes(x),
  indiv_missingness = function(x, g, xc) indiv_missingness(x),
  loci_missingness = function(x, g, xc) loci_missingness(x),
  loci_missingness_grouped = function(x, g, xc) loci_missingness(g),
  loci_alt_freq = function(x, g, xc) loci_alt_freq(x),
  loci_alt_freq_grouped = function(x, g, xc) loci_alt_freq(g),
  loci_maf = function(x, g, xc) loci_maf(x),
  loci_maf_grouped = function(x, g, xc) loci_maf(g),
  loci_hwe = function(x, g, xc) loci_hwe(x),
  loci_hwe_grouped = function(x, g, xc) loci_hwe(g),
  loci_pi = function(x, g, xc) loci_pi(x),
  loci_pi_grouped = function(x, g, xc) loci_pi(g),
  loci_ld_clump = function(x, g, xc) loci_ld_clump(xc),
  select_loci_if = function(x, g, xc) {
    select_loci_if(x, loci_maf(genotypes) > 0.05)
  },
  # individual statistics
  indiv_het_obs = function(x, g, xc) indiv_het_obs(x),
  indiv_inbreeding = function(x, g, xc) indiv_inbreeding(x),
  windows_indiv_roh = function(x, g, xc) {
    windows_indiv_roh(x, window_snp = 10, min_snp = 5, min_kb = 1)
  },
  # population statistics
  pop_het_obs = function(x, g, xc) pop_het_obs(g),
  pop_het_exp = function(x, g, xc) pop_het_exp(g),
  pop_gene_div = function(x, g, xc) pop_gene_div(g),
  pop_fis_Nei87 = function(x, g, xc) pop_fis(g, method = "Nei87"),
  pop_fis_WG17 = function(x, g, xc) pop_fis(g, method = "WG17"),
  pop_fst = function(x, g, xc) pop_fst(g),
  pop_global_stats = function(x, g, xc) pop_global_stats(g),
  pop_tajimas_d = function(x, g, xc) pop_tajimas_d(g),
  windows_pop_tajimas_d = function(x, g, xc) {
    windows_pop_tajimas_d(g, window_size = 50, step_size = 50)
  },
  pop_sfs = function(x, g, xc) pop_sfs(g),
  pairwise_pop_sfs = function(x, g, xc) pairwise_pop_sfs(g),
  pairwise_pop_fst_Hudson = function(x, g, xc) {
    pairwise_pop_fst(g, method = "Hudson")
  },
  pairwise_pop_fst_Nei87 = function(x, g, xc) pairwise_pop_fst(g, method = "Nei87"),
  pairwise_pop_fst_WC84 = function(x, g, xc) pairwise_pop_fst(g, method = "WC84"),
  windows_pairwise_pop_fst = function(x, g, xc) {
    windows_pairwise_pop_fst(g, window_size = 50, step_size = 50)
  },
  nwise_pop_pbs = function(x, g, xc) nwise_pop_pbs(g),
  # pairwise individual statistics
  pairwise_ibs = function(x, g, xc) pairwise_ibs(x),
  pairwise_allele_sharing = function(x, g, xc) pairwise_allele_sharing(x),
  pairwise_grm = function(x, g, xc) pairwise_grm(x),
  pairwise_king = function(x, g, xc) pairwise_king(x),
  # qc reports
  qc_report_loci = function(x, g, xc) qc_report_loci(x),
  qc_report_loci_grouped = function(x, g, xc) qc_report_loci(g),
  qc_report_indiv = function(x, g, xc) qc_report_indiv(x),
  # imputation and multivariate analyses
  gt_impute_simple = function(x, g, xc) gt_impute_simple(x),
  gt_pca_partialSVD = function(x, g, xc) gt_pca_partialSVD(xc),
  gt_pca_randomSVD = function(x, g, xc) gt_pca_randomSVD(xc),
  gt_pca_autoSVD = function(x, g, xc) {
    gt_pca_autoSVD(xc, k = 3, roll_size = 5)
  },
  gt_dapc = function(x, g, xc) {
    pca <- gt_pca_partialSVD(xc)
    gt_dapc(pca, pop = factor(xc$population), n_pca = 5)
  },
  # manipulation
  gt_pseudohaploid = function(x, g, xc) gt_pseudohaploid(x),
  rbind = function(x, g, xc) {
    rbind(x[1:18, ], x[19:36, ], backingfile = tempfile("rb_"), quiet = TRUE)
  },
  gt_update_backingfile = function(x, g, xc) {
    gt_update_backingfile(x[1:10, ], backingfile = tempfile("ub_"),
                          quiet = TRUE)
  },
  gt_save_load = function(x, g, xc) {
    f <- gt_save(x, file_name = tempfile(fileext = ".gt"), quiet = TRUE)
    gt_load(f[1])
  },
  filter_refresh = function(x, g, xc) show_ploidy(x[indiv_ploidy(x) ==
                                                  max(indiv_ploidy(x)), ]),
  # export
  gt_as_vcf = function(x, g, xc) gt_as_vcf(x, file = tempfile(fileext = ".vcf")),
  gt_as_plink = function(x, g, xc) gt_as_plink(x, file = tempfile(fileext = ".bed")),
  gt_as_genind = function(x, g, xc) gt_as_genind(x),
  gt_as_genlight = function(x, g, xc) gt_as_genlight(x),
  gt_as_hierfstat = function(x, g, xc) gt_as_hierfstat(x),
  gt_as_geno_lea = function(x, g, xc) gt_as_geno_lea(x),
  gt_extract_f2 = function(x, g, xc) {
    gt_extract_f2(g, outdir = tempfile("f2_"), quiet = TRUE)
  }
)

run_call <- function(f, x, g, xc) {
  warns <- character()
  res <- withCallingHandlers(
    tryCatch(
      {
        f(x, g, xc)
        list(status = "ok", message = "")
      },
      error = function(e) {
        list(status = "error", message = conditionMessage(e))
      }
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # a suggested package that is not installed says nothing about ploidy
  if (grepl("first install package", res$message)) {
    res$status <- "no_pkg"
  }
  if (res$status == "ok" && length(warns) > 0) {
    res <- list(status = "warning", message = warns[1])
  }
  res$message <- substr(gsub("\\s+", " ", res$message), 1, 160)
  res
}

rows <- list()
for (d in names(datasets)) {
  x <- datasets[[d]]
  g <- group_by(x, population)
  for (fn in names(calls)) {
    r <- suppressMessages(run_call(calls[[fn]], x, g, datasets_complete[[d]]))
    rows[[length(rows) + 1]] <- data.frame(
      call = fn, dataset = d, status = r$status, message = r$message
    )
  }
}
long <- do.call(rbind, rows)
write.table(long, sub("\\.tsv$", "_long.tsv", out_file),
  sep = "\t", quote = FALSE, row.names = FALSE
)
wide <- tidyr::pivot_wider(long[, 1:3], names_from = "dataset",
                           values_from = "status")
write.table(wide, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(wide))
