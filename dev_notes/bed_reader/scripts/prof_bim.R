# Time the steps of loci_from_bim() after fread
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(tidypopgen, lib.loc = args[1]))
ns <- asNamespace("tidypopgen")
bim <- args[2]
tm <- function(label, expr) {
  t0 <- proc.time(); res <- force(expr); t1 <- proc.time() - t0
  cat(sprintf("%-34s elapsed %6.2f\n", label, t1[["elapsed"]])); invisible(res)
}
b <- tm("fread", ns$read_loci_file(bim, character_cols = c(2, 5, 6)))
names(b)[1:6] <- c("chromosome", "marker.ID", "genetic.dist", "physical.pos", "allele1", "allele2")
loci <- tm("tibble + as.factor etc", tibble::tibble(
  name = as.character(b$marker.ID), chromosome = as.factor(b$chromosome),
  position = as.integer(b$physical.pos), genetic_dist = as.numeric(b$genetic.dist),
  allele_ref = as.character(b$allele2), allele_alt = as.character(b$allele1)))
tm("  as.factor(chromosome) alone", as.factor(b$chromosome))
tm("check_allele_alphabet", ns$check_allele_alphabet(loci))
loci <- tm("harmonise_missing_values", ns$harmonise_missing_values(loci))
cg <- tm("chromosome_group_codes", ns$chromosome_group_codes(loci$chromosome))
tm("unsorted_by_chromosome_cpp", ns$unsorted_by_chromosome_cpp(loci$position, cg$codes, seq_len(cg$n)))
tm("any_duplicated_names", ns$any_duplicated_names(loci$name))
tm("validate_loci (all)", ns$validate_loci(loci, check_alphabet = TRUE, check_duplicates = TRUE, harmonise_loci = TRUE))
