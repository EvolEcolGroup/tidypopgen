# The tidypopgen side of the VCFArrow comparison
# (https://github.com/legalLab/VCFArrow/blob/main/comparison/tool_comparison.md).
# VCFArrow did not publish its scripts, so each step is our reading of the
# call named in their document ("What each tool does for the same function");
# `vcfarrow` gives the row of their tables, and `extra` marks steps that are
# not in their tables (added to look for weak spots).
#
# Each step is function(gt, ctx) and returns a short character summary
# (dimensions, counts, file sizes), so that the versions can be checked to
# give the same answer. Everything runs from the unfiltered data, as in the
# comparison. ctx: threads, out_dir, strata, vcf, work.

suppressPackageStartupMessages({
  library(dplyr)
})

dims <- function(x) {
  sprintf("%d indiv x %d loci", nrow(x), count_loci(x))
}
fsize <- function(files) {
  files <- files[file.exists(files)]
  sprintf("%s: %.1f MB", basename(files), file.size(files) / 1e6) |>
    paste(collapse = "; ")
}
num_summary <- function(v) {
  sprintf("n=%d, NA=%d, mean=%.6g", length(v), sum(is.na(v)),
    mean(v, na.rm = TRUE))
}
# remove a gen_tibble's output files (big at full scale)
rm_out <- function(ctx, pattern) {
  unlink(list.files(ctx$out_dir, pattern = pattern, full.names = TRUE))
}

steps <- list(
  # read_vcf and set_groups run in the same process (see run_step.R)
  read_vcf = list(vcfarrow = "read VCF", fun = function(gt, ctx) {
    gen_tibble(ctx$vcf,
      backingfile = file.path(ctx$work, "gt"),
      parser = "cpp", n_cores = ctx$threads, quiet = TRUE
    )
  }),
  set_groups = list(vcfarrow = "set groups", fun = function(gt, ctx) {
    strata <- utils::read.delim(ctx$strata)
    gt %>% mutate(population = strata$pop[match(id, strata$id)])
  }),
  extract_samples = list(vcfarrow = "extract samples", fun = function(gt, ctx) {
    dims(gt[1:10, ])
  }),
  extract_groups = list(vcfarrow = "extract groups", fun = function(gt, ctx) {
    dims(filter(gt, population %in% c("GS", "BS", "WA")))
  }),
  subsample_loci = list(vcfarrow = "subsample loci (1000)", fun = function(gt, ctx) {
    set.seed(1)
    idx <- sort(sample.int(count_loci(gt), 1000))
    dims(select_loci(gt, .sel_arg = all_of(idx)))
  }),
  filter_maf = list(vcfarrow = "filter MAF (0.05)", fun = function(gt, ctx) {
    nc <- ctx$threads
    dims(select_loci_if(gt, coalesce(loci_maf(genotypes, n_cores = nc) >= 0.05, FALSE)))
  }),
  filter_locus_missing = list(vcfarrow = "filter locus missing (0.2)", fun = function(gt, ctx) {
    nc <- ctx$threads
    dims(select_loci_if(gt, loci_missingness(genotypes, n_cores = nc) <= 0.2))
  }),
  filter_sample_missing = list(vcfarrow = "filter sample missing (0.8)", fun = function(gt, ctx) {
    dims(filter(gt, indiv_missingness(genotypes) <= 0.8))
  }),
  filter_invariant = list(vcfarrow = "filter invariant", fun = function(gt, ctx) {
    nc <- ctx$threads
    dims(select_loci_if(gt, coalesce(loci_maf(genotypes, n_cores = nc) > 0, FALSE)))
  }),
  merge = list(vcfarrow = "merge two VCFs", fun = function(gt, ctx) {
    # two sample sets of the same data (ingroups and outgroup), merged into
    # a new backing file
    ingroup <- filter(gt, population != "OG")
    outgroup <- filter(gt, population == "OG")
    merged <- rbind(ingroup, outgroup,
      backingfile = file.path(ctx$out_dir, "merged"), quiet = TRUE
    )
    res <- dims(merged)
    rm_out(ctx, "^merged")
    res
  }),
  sample_stats = list(vcfarrow = "sample statistics", fun = function(gt, ctx) {
    qc <- qc_report_indiv(gt)
    sprintf("het_obs %s; missingness %s", num_summary(qc$het_obs),
      num_summary(qc$missingness))
  }),
  theta_pi = list(vcfarrow = "theta / pi", fun = function(gt, ctx) {
    num_summary(loci_pi(gt, n_cores = ctx$threads))
  }),
  missing_plot = list(vcfarrow = "missing data plot + table", fun = function(gt, ctx) {
    miss <- tibble::tibble(id = gt$id, population = gt$population,
      missingness = indiv_missingness(gt$genotypes))
    utils::write.csv(miss, file.path(ctx$out_dir, "missing.csv"),
      row.names = FALSE)
    p <- ggplot2::ggplot(miss, ggplot2::aes(id, missingness, fill = population)) +
      ggplot2::geom_col()
    ggplot2::ggsave(file.path(ctx$out_dir, "missing.png"), p,
      width = 8, height = 4)
    num_summary(miss$missingness)
  }),
  write_vcf = list(vcfarrow = "write VCF", fun = function(gt, ctx) {
    f <- file.path(ctx$out_dir, "out.vcf")
    gt_as_vcf(gt, file = f, overwrite = TRUE)
    res <- fsize(f)
    rm_out(ctx, "^out\\.vcf")
    res
  }),
  to_genlight = list(vcfarrow = "-> genlight", fun = function(gt, ctx) {
    gl <- gt_as_genlight(gt)
    sprintf("genlight %d x %d", adegenet::nInd(gl), adegenet::nLoc(gl))
  }),
  to_ped = list(vcfarrow = "-> PLINK ped", fun = function(gt, ctx) {
    f <- file.path(ctx$out_dir, "out_ped")
    gt_as_plink(gt, file = f, type = "ped", overwrite = TRUE)
    res <- fsize(paste0(f, c(".ped", ".map")))
    rm_out(ctx, "^out_ped")
    res
  }),
  to_bed = list(vcfarrow = "-> PLINK bed", fun = function(gt, ctx) {
    f <- file.path(ctx$out_dir, "out_bed")
    gt_as_plink(gt, file = f, type = "bed", overwrite = TRUE)
    res <- fsize(paste0(f, c(".bed", ".bim", ".fam")))
    rm_out(ctx, "^out_bed")
    res
  }),

  # extra steps -----------------------------------------------------------
  # the same subsample, with the index vector given directly (no all_of())
  subsample_loci_bare = list(vcfarrow = NA, extra = TRUE, fun = function(gt, ctx) {
    set.seed(1)
    idx <- sort(sample.int(count_loci(gt), 1000))
    dims(suppressWarnings(select_loci(gt, .sel_arg = idx)))
  }),
  # pi per group, closer to what vcf_theta() reports
  theta_pi_grouped = list(vcfarrow = NA, extra = TRUE, fun = function(gt, ctx) {
    pi <- gt %>%
      group_by(population) %>%
      loci_pi(n_cores = ctx$threads, type = "matrix")
    sprintf("%d x %d, mean=%.6g", nrow(pi), ncol(pi), mean(pi, na.rm = TRUE))
  }),
  # the filter chain of the VCFArrow vignette, then written to a new backing
  # file, as VCFArrow filters always produce a complete new object
  filter_chain_update = list(vcfarrow = NA, extra = TRUE, fun = function(gt, ctx) {
    nc <- ctx$threads
    x <- gt %>%
      filter(population != "OG") %>%
      select_loci_if(coalesce(loci_maf(genotypes, n_cores = nc) >= 0.03, FALSE)) %>%
      select_loci_if(loci_missingness(genotypes, n_cores = nc) <= 0.2) %>%
      filter(indiv_missingness(genotypes) <= 0.3)
    x <- gt_update_backingfile(x,
      backingfile = file.path(ctx$out_dir, "filtered"), quiet = TRUE
    )
    res <- dims(x)
    rm_out(ctx, "^filtered")
    res
  })
)

# read with integer loci names (1.0.0 and later); run as its own step, as it
# writes a second backing file
steps$read_vcf_int <- list(vcfarrow = NA, extra = TRUE, needs = "names_as_int",
  fun = function(gt, ctx) {
    x <- gen_tibble(ctx$vcf,
      backingfile = file.path(ctx$out_dir, "gt_int"),
      parser = "cpp", n_cores = ctx$threads, names_as_int = TRUE, quiet = TRUE
    )
    res <- dims(x)
    rm_out(ctx, "^gt_int")
    res
  })
