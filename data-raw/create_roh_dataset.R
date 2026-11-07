# Create a small simulated dataset with known runs of homozygosity, to test
# and illustrate windows_indiv_roh() without external data.
#
# Writes, in inst/extdata/roh:
# - roh_sim.bed/.bim/.fam: 30 diploid individuals in three populations, and
#   3,000 loci on three chromosomes (1,000 loci each, over 5, 4 and 3 Mb, with
#   genetic distances of 1 cM per Mb);
# - roh_sim_truth.csv: the simulated runs (one row per run), with the
#   individual, the chromosome, the rows of the first and last locus of the
#   run in the loci table, their positions, the number of loci, and the
#   number of heterozygous and missing genotypes in the run.
#
# Each individual is made of two haplotypes, drawn independently from the
# allele frequencies of each locus (uniform between 0.3 and 0.7, so that
# windows of unrelated haplotypes are rarely homozygous). In a run of
# homozygosity (an autozygous segment), the second haplotype is a copy of the
# first. The populations differ in inbreeding:
# - "outbred": no runs;
# - "moderate": one run of 60 to 150 loci per individual;
# - "inbred": three runs of 150 to 350 loci per individual, except for the
#   first three individuals, which carry runs designed to test edge cases:
#   - inbred_01: a run at the end of chromosome 1 (to its last locus), and
#     one at the start of chromosome 2 (from its first locus); they must be
#     two runs, one per chromosome;
#   - inbred_02: a run on chromosome 3 with one heterozygous genotype (a
#     genotyping error) and two missing genotypes, more than 15 loci apart,
#     so that with the default settings it is a single run;
#   - inbred_03: a run on chromosome 2 that spans a gap of 1.3 Mb without
#     loci, larger than the default max_gap of windows_indiv_roh(), so that
#     it is found as two runs, one on each side of the gap.
# About 0.2% of all genotypes are missing at random.
#
# Run from the root of the repository.

set.seed(1234)
out_dir <- "inst/extdata/roh"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# loci ------------------------------------------------------------------------
n_per_chrom <- 1000
chrom_length <- c(5e6, 4e6, 3e6)
# chromosome 2 has no loci between 1.5 and 2.8 Mb
gap <- c(1.5e6, 2.8e6)
positions <- lapply(seq_along(chrom_length), function(chrom) {
  candidates <- seq_len(chrom_length[chrom])
  if (chrom == 2) {
    candidates <- candidates[candidates <= gap[1] | candidates > gap[2]]
  }
  sort(sample(candidates, n_per_chrom))
})
n_loci <- n_per_chrom * length(chrom_length)
alleles <- t(replicate(n_loci, sample(c("A", "C", "G", "T"), 2)))
loci <- data.frame(
  name = paste0("snp", seq_len(n_loci)),
  chromosome = rep(seq_along(chrom_length), each = n_per_chrom),
  position = as.integer(unlist(positions)),
  genetic_dist = unlist(positions) / 1e6,
  allele_ref = alleles[, 1],
  allele_alt = alleles[, 2]
)
# the row of a locus in the loci table, from its chromosome and its index on
# that chromosome
locus_row <- function(chrom, i) (chrom - 1) * n_per_chrom + i
# the last locus before the gap on chromosome 2
gap_index <- max(which(positions[[2]] <= gap[1]))

# individuals -----------------------------------------------------------------
n_per_pop <- 10
pops <- c("outbred", "moderate", "inbred")
indiv_meta <- data.frame(
  id = paste0(rep(pops, each = n_per_pop), "_",
              sprintf("%02d", seq_len(n_per_pop))),
  population = rep(pops, each = n_per_pop)
)

# runs ------------------------------------------------------------------------
# a run is given by its chromosome and the indices of its first and last loci
# on that chromosome
run <- function(id, chrom, from, to) {
  data.frame(id = id, chrom = chrom, from = from, to = to)
}
# random runs that do not overlap (with a margin of 50 loci) each other, the
# ends of the chromosomes, or the gap on chromosome 2
random_runs <- function(id, n, min_len, max_len) {
  runs <- run(character(0), integer(0), integer(0), integer(0))
  while (nrow(runs) < n) {
    chrom <- sample(seq_along(chrom_length), 1)
    len <- sample(min_len:max_len, 1)
    from <- sample(50:(n_per_chrom - len - 50), 1)
    to <- from + len - 1
    same_chrom <- runs$chrom == chrom
    overlaps <- any(from - 50 <= runs$to[same_chrom] &
                      to + 50 >= runs$from[same_chrom])
    on_gap <- chrom == 2 && from - 50 <= gap_index + 1 &&
      to + 50 >= gap_index
    if (!overlaps && !on_gap) runs <- rbind(runs, run(id, chrom, from, to))
  }
  runs
}
ids <- split(indiv_meta$id, indiv_meta$population)
runs <- rbind(
  do.call(rbind, lapply(ids$moderate, random_runs, 1, 60, 150)),
  run("inbred_01", 1, 801, n_per_chrom),
  run("inbred_01", 2, 1, 150),
  run("inbred_02", 3, 300, 600),
  run("inbred_03", 2, gap_index - 119, gap_index + 120),
  do.call(rbind, lapply(ids$inbred[-(1:3)], random_runs, 3, 150, 350))
)
runs$first_locus <- locus_row(runs$chrom, runs$from)
runs$last_locus <- locus_row(runs$chrom, runs$to)

# genotypes -------------------------------------------------------------------
freq <- stats::runif(n_loci, 0.3, 0.7)
haplotype <- function() stats::rbinom(n_loci, 1, freq)
genotypes <- t(vapply(indiv_meta$id, function(id) {
  hap_1 <- haplotype()
  hap_2 <- haplotype()
  for (i in which(runs$id == id)) {
    in_run <- runs$first_locus[i]:runs$last_locus[i]
    hap_2[in_run] <- hap_1[in_run]
  }
  hap_1 + hap_2
}, numeric(n_loci)))
dimnames(genotypes) <- NULL
# missing genotypes at random
genotypes[sample(length(genotypes), round(0.002 * length(genotypes)))] <- NA
# the heterozygous and missing genotypes in the run of inbred_02
inbred_02 <- which(indiv_meta$id == "inbred_02")
genotypes[inbred_02, locus_row(3, 450)] <- 1
genotypes[inbred_02, locus_row(3, c(400, 500))] <- NA

# write the files -------------------------------------------------------------
sim_gt <- tidypopgen::gen_tibble(
  genotypes,
  indiv_meta = indiv_meta,
  loci = loci,
  backingfile = tempfile(),
  quiet = TRUE
)
tidypopgen::gt_as_plink(sim_gt, file = file.path(out_dir, "roh_sim"))

truth <- data.frame(
  id = runs$id,
  population = indiv_meta$population[match(runs$id, indiv_meta$id)],
  chromosome = runs$chrom,
  first_locus = runs$first_locus,
  last_locus = runs$last_locus,
  start = loci$position[runs$first_locus],
  end = loci$position[runs$last_locus],
  n_loci = runs$last_locus - runs$first_locus + 1,
  n_het = mapply(function(id, from, to) {
    sum(genotypes[indiv_meta$id == id, from:to] == 1, na.rm = TRUE)
  }, runs$id, runs$first_locus, runs$last_locus),
  n_missing = mapply(function(id, from, to) {
    sum(is.na(genotypes[indiv_meta$id == id, from:to]))
  }, runs$id, runs$first_locus, runs$last_locus)
)
truth <- truth[order(match(truth$id, indiv_meta$id), truth$first_locus), ]
utils::write.csv(truth, file.path(out_dir, "roh_sim_truth.csv"),
  row.names = FALSE
)
