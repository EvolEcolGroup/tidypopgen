# Rscript edge_check.R <lib> <out_rds>: merges small VCFs with swapped loci,
# haploid genotypes in diploid individuals, missing genotypes, phased and
# triploid files, and saves the genotypes, to compare versions
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(tidypopgen, lib.loc = args[1]))
d <- file.path(Sys.getenv("TPG_WORK", "/media/andrea/Elements1/tmp/claude/tidypopgen"), "merge_vcf/edge")
hdr <- function(s) c("##fileformat=VCFv4.2", paste(c("#CHROM", "POS", "ID",
  "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", s), collapse = "\t"))
line <- function(chr, pos, ref, alt, gts) {
  paste(c(chr, pos, ".", ref, alt, ".", "PASS", ".", "GT", gts), collapse = "\t")
}
write_vcf <- function(name, samples, lines) {
  f <- file.path(d, name)
  writeLines(c(hdr(samples), lines), f)
  f
}
a <- write_vcf("a.vcf", c("a1", "a2"), c(
  line(1, 100, "A", "C", c("0/1", "1/1")),
  line(1, 200, "G", "T", c("0/0", "./.")),
  line("X", 300, "A", "G", c("0", "1")),
  line("X", 400, "C", "T", c("0/1", "1"))
))
b <- write_vcf("b.vcf", c("b1", "b2"), c(
  line(1, 100, "C", "A", c("0/1", "0/0")), # swapped
  line(1, 200, "G", "T", c("1/1", "0/1")),
  line("X", 300, "G", "A", c("1", "0")), # swapped, haploid
  line("X", 400, "C", "T", c("0", "./."))
))
c3 <- write_vcf("c.vcf", c("c1"), c(
  line(1, 100, "A", "C", c("0/1/1")),
  line(1, 200, "T", "G", c("1/1/0")), # swapped, triploid
  line("X", 300, "A", "G", c("0/0/1")),
  line("X", 400, "C", "T", c("././."))
))
p1 <- write_vcf("p1.vcf", c("p1"), c(
  line(1, 100, "A", "C", c("0|1")),
  line(1, 200, "G", "T", c("1|0")),
  line("X", 300, "A", "G", c(".|.")),
  line("X", 400, "C", "T", c("1|1"))
))
p2 <- write_vcf("p2.vcf", c("p2"), c(
  line(1, 100, "A", "C", c("1|1")),
  line(1, 200, "G", "T", c("0|1")),
  line("X", 400, "C", "T", c("0|0"))
))
res <- list()
get_gt <- function(x) {
  list(g = show_genotypes(x), loci = show_loci(x)[, -1], ploidy = indiv_ploidy(x))
}
res$ab <- get_gt(gt_merge_vcfs(c(a, b), ref_mismatch = "swap",
  backingfile = tempfile(), quiet = TRUE))
res$abc <- get_gt(gt_merge_vcfs(c(a, b, c3), ref_mismatch = "swap",
  backingfile = tempfile(), quiet = TRUE))
res$abc_union <- get_gt(gt_merge_vcfs(c(c3, a, b), ref_mismatch = "swap",
  variants = "union", backingfile = tempfile(), quiet = TRUE))
res$p <- get_gt(gt_merge_vcfs(c(p1, p2), variants = "union", phased = TRUE,
  backingfile = tempfile(), quiet = TRUE))
res$p_ref <- get_gt(gt_merge_vcfs(c(p1, p2), variants = "union",
  missing_to_ref = TRUE, backingfile = tempfile(), quiet = TRUE))
saveRDS(res, args[2])
print(res$abc$g)
print(res$p$g)
