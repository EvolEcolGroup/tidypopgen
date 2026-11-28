# Test fixture for indiv_sex_check() and the coding of hemizygous genotypes:
# a small simulated dataset with autosomes, X, Y and MT, and the outputs of
# PLINK 1.9 on it (--check-sex, with and without ycount), and PLINK's .bed
# made from a VCF of the same genotypes in which males are haploid on X, Y
# and MT (to check that gen_tibble() reads haploid genotypes as PLINK).
#
# 18 individuals: 8 females and 8 males (founders), a daughter of the first
# female and male (not a founder, so not used for allele frequencies), and
# a female of unknown sex. Two individuals are misrecorded: a female
# recorded as male, and a male recorded as female. One male has a few
# heterozygous genotypes on X, one female a few genotypes on Y. Loci: 20 on
# chromosome 1, 200 on X (one monomorphic), 10 on Y, 5 on MT; 3% of the
# genotypes are missing.
#
# Run from the package root, with PLINK 1.9 on the PATH (or in PLINK):
#   Rscript data-raw/validate/create_plink_sexcheck_test_data.R
devtools::load_all(quiet = TRUE)
out_dir <- "tests/testthat/testdata/plink_sexcheck"
dir.create(out_dir, showWarnings = FALSE)
plink <- Sys.getenv("PLINK", "plink")

set.seed(123)
# the genetic sex of each individual (1 male, 2 female), and its recorded sex
# (PLINK's codes, 0 unknown)
genetic_sex <- c(rep(2, 8), rep(1, 8), 2, 2)
recorded_sex <- genetic_sex
recorded_sex[3] <- 1 # a female recorded as male
recorded_sex[12] <- 2 # a male recorded as female
recorded_sex[18] <- 0 # a female of unknown sex
n_ind <- length(genetic_sex)
male <- genetic_sex == 1
n_auto <- 20
n_x <- 200
n_y <- 10
n_mt <- 5
# diploid genotypes, or haploid ones stored as homozygotes (0 or 2)
diploid <- function(p, n) t(replicate(n, stats::rbinom(length(p), 2, p)))
haploid <- function(p, n) 2 * t(replicate(n, stats::rbinom(length(p), 1, p)))
p_auto <- stats::runif(n_auto, 0.1, 0.9)
p_x <- stats::runif(n_x, 0.1, 0.9)
p_y <- stats::runif(n_y, 0.2, 0.8)
p_mt <- stats::runif(n_mt, 0.2, 0.8)
auto <- diploid(p_auto, n_ind)
x <- matrix(NA, n_ind, n_x)
x[!male, ] <- diploid(p_x, sum(!male))
x[male, ] <- haploid(p_x, sum(male))
x[, 1] <- 0 # a monomorphic locus
# the daughter (17) of individuals 1 and 9 takes her father's X
x[17, ] <- x[9, ]
x[17, x[1, ] == 2] <- x[17, x[1, ] == 2] / 2 + 1
y <- matrix(NA, n_ind, n_y)
y[male, ] <- haploid(p_y, sum(male))
y[5, 1:3] <- 0 # a female with a few genotypes on Y
mt <- haploid(p_mt, n_ind)
# a male with a few heterozygous genotypes on X
x[10, 2:4] <- 1
geno <- cbind(auto, x, y, mt)
geno[matrix(stats::runif(length(geno)), nrow(geno)) < 0.03] <- NA

chromosome <- rep(c("1", "X", "Y", "MT"), c(n_auto, n_x, n_y, n_mt))
loci <- data.frame(
  name = paste0("snp", seq_len(ncol(geno))), chromosome = chromosome,
  position = unlist(lapply(c(n_auto, n_x, n_y, n_mt), seq_len)) * 1000L,
  genetic_dist = 0, allele_ref = "A", allele_alt = "G"
)
indiv_meta <- data.frame(
  id = paste0("ind", seq_len(n_ind)),
  population = "pop1",
  sex = c("male", "female")[ifelse(recorded_sex == 0, NA, recorded_sex)]
)
gt <- gen_tibble(geno,
  indiv_meta = indiv_meta, loci = loci,
  backingfile = tempfile(), quiet = TRUE
)
bfile <- file.path(out_dir, "sexcheck")
gt_as_plink(gt, file = paste0(bfile, ".bed"), type = "bed")
# the parents of the daughter (gt_as_plink() writes 0 for unknown parents)
fam <- utils::read.table(paste0(bfile, ".fam"))
fam[17, 3:4] <- c("ind9", "ind1")
utils::write.table(fam, paste0(bfile, ".fam"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

# PLINK's sex checks: default thresholds, and ycount with males needing at
# least 5 genotypes on Y
system2(plink, c("--bfile", bfile, "--check-sex", "--out", bfile))
system2(plink, c(
  "--bfile", bfile, "--check-sex", "ycount", "0.2", "0.8", "0", "5",
  "--out", paste0(bfile, "_ycount")
))

# a VCF of the same genotypes, with males haploid on X, Y and MT, and
# everyone haploid on MT (written by hand, so that it does not depend on
# gt_as_vcf()); females have no genotypes on Y, written as "."
vcf_gt <- function(g, haploid) {
  ifelse(is.na(g), if (haploid) "." else "./.",
    if (haploid) {
      ifelse(g == 0, "0", ifelse(g == 2, "1", "0/1"))
    } else {
      c("0/0", "0/1", "1/1")[g + 1]
    }
  )
}
vcf_lines <- vapply(seq_len(ncol(geno)), function(j) {
  hap <- switch(chromosome[j],
    "1" = rep(FALSE, n_ind),
    "X" = male,
    "Y" = rep(TRUE, n_ind),
    "MT" = rep(TRUE, n_ind)
  )
  gts <- vapply(seq_len(n_ind), function(i) {
    vcf_gt(geno[i, j], hap[i])
  }, character(1))
  paste(c(
    chromosome[j], loci$position[j], loci$name[j], "A", "G", ".", "PASS",
    ".", "GT", gts
  ), collapse = "\t")
}, character(1))
vcf_file <- paste0(bfile, ".vcf")
writeLines(c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste(c(
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
    indiv_meta$id
  ), collapse = "\t"),
  vcf_lines
), vcf_file)
# PLINK's .bed of the VCF (keeping the alleles in the order of the VCF)
system2(plink, c(
  "--vcf", vcf_file, "--keep-allele-order", "--make-bed",
  "--out", paste0(bfile, "_vcf")
))

writeLines(system2(plink, "--version", stdout = TRUE),
  file.path(out_dir, "plink_version.txt")
)
# keep only the files the tests use
unlink(paste0(bfile, c(".log", ".nosex", ".hh")))
unlink(paste0(bfile, "_ycount", c(".log", ".nosex", ".hh")))
unlink(paste0(bfile, "_vcf", c(".log", ".nosex", ".hh")))
