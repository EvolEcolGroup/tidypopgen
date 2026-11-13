# A synthetic PLINK .ped file (with its .map) to benchmark gen_tibble() on
# .ped files, written by gt_as_plink() from a random gen_tibble.
#   Rscript make_data.R <lib> <out_dir> [n_indiv] [n_loci]
# (lib is the library holding the tidypopgen used to write the file)
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
out_dir <- args[2]
n_indiv <- if (length(args) >= 3) as.integer(args[3]) else 2000L
n_loci <- if (length(args) >= 4) as.integer(args[4]) else 50000L
suppressPackageStartupMessages(library(tidypopgen, lib.loc = lib))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
# dosages with locus-specific allele frequencies, and 2% missing genotypes
freq <- stats::runif(n_loci, 0.05, 0.95)
geno <- matrix(
  stats::rbinom(n_indiv * n_loci, 2, rep(freq, each = n_indiv)),
  nrow = n_indiv
)
geno[sample(length(geno), round(0.02 * length(geno)))] <- NA
bases <- c("A", "C", "G", "T")
ref <- sample(bases, n_loci, replace = TRUE)
alt <- vapply(ref, function(r) sample(setdiff(bases, r), 1), character(1))
loci <- data.frame(
  name = paste0("rs", seq_len(n_loci)),
  chromosome = rep(1:20, length.out = n_loci) |> sort(),
  position = seq_len(n_loci) * 100L,
  genetic_dist = 0,
  allele_ref = ref,
  allele_alt = unname(alt)
)
indiv_meta <- data.frame(
  id = paste0("ind", seq_len(n_indiv)),
  population = paste0("pop", rep(1:10, length.out = n_indiv))
)
gt <- gen_tibble(geno,
  indiv_meta = indiv_meta, loci = loci,
  backingfile = file.path(out_dir, "synthetic_src"), quiet = TRUE
)
rm(geno)
ped <- gt_as_plink(gt, file = file.path(out_dir, "synthetic.ped"), type = "ped")
cat("written", ped, round(file.size(ped) / 2^20), "MB\n")
