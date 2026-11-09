# write a small VCF (with GT only) for the tests: loci is a data.frame with
# chrom, pos, id, ref and alt, and genotypes a character matrix (one row per
# locus, one column per sample)
write_merge_vcf <- function(samples, loci, genotypes,
                            file = tempfile(fileext = ".vcf")) {
  header <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", samples
    ), collapse = "\t")
  )
  body <- paste(loci$chrom, loci$pos, loci$id, loci$ref, loci$alt, ".",
    "PASS", ".", "GT", apply(genotypes, 1, paste, collapse = "\t"),
    sep = "\t"
  )
  writeLines(c(header, body), file)
  file
}

# two files with different samples:
# - 1:100 has a different ALT (a multiallelic site, dropped by default)
# - 1:200 is the same in both files
# - 1:300 is monomorphic in file 2 (ALT "."), so it takes the ALT of file 1
# - 2:50 is only in file 1, 1:400 only in file 2
vcf_1 <- write_merge_vcf(
  samples = c("a", "b"),
  loci = data.frame(
    chrom = c("1", "1", "1", "2"), pos = c(100, 200, 300, 50),
    id = c("s1", "s2", "s3", "s4"), ref = c("A", "G", "C", "T"),
    alt = c("C", "T", "G", "A")
  ),
  genotypes = rbind(
    c("0/1", "1/1"), c("0/0", "0/1"), c("1/1", "./."), c("0/1", "0/0")
  )
)
vcf_2 <- write_merge_vcf(
  samples = "c",
  loci = data.frame(
    chrom = c("1", "1", "1", "1"), pos = c(100, 200, 300, 400),
    id = c("s1", "s2", "s3", "s5"), ref = c("A", "G", "C", "A"),
    alt = c("G", "T", ".", "T")
  ),
  genotypes = rbind("1/1", "0/1", "0/0", "0/1")
)

test_that("gt_merge_vcfs keeps the intersection or the union of the sites", {
  # intersection (the default): sites in both files
  merged <- gt_merge_vcfs(c(vcf_1, vcf_2),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(merged$id, c("a", "b", "c"))
  expect_identical(show_loci(merged)$name, c("s2", "s3"))
  expect_identical(show_loci(merged)$allele_alt, c("T", "G"))
  expect_equal(
    show_genotypes(merged),
    cbind(c(0, 1, 1), c(2, NA, 0))
  )
  # the merged gen_tibble is saved, and can be loaded back
  reloaded <- gt_load(bigstatsr::sub_bk(gt_get_file_names(merged)[2], ".gt"))
  expect_equal(show_genotypes(reloaded), show_genotypes(merged))

  # union: all sites but the multiallelic one, sorted by chromosome (in order
  # of appearance) and position, with missing genotypes for the samples of
  # the file without the site
  merged <- gt_merge_vcfs(c(vcf_1, vcf_2),
    variants = "union",
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(show_loci(merged)$name, c("s2", "s3", "s5", "s4"))
  expect_identical(
    as.character(show_loci(merged)$chromosome),
    c("1", "1", "1", "2")
  )
  expect_equal(
    show_genotypes(merged),
    cbind(c(0, 1, 1), c(2, NA, 0), c(NA, NA, 1), c(1, 0, NA))
  )
  # or homozygous for the reference allele, as with bcftools merge -0
  merged <- gt_merge_vcfs(c(vcf_1, vcf_2),
    variants = "union", missing_to_ref = TRUE,
    backingfile = tempfile(), quiet = TRUE
  )
  expect_equal(
    show_genotypes(merged),
    cbind(c(0, 1, 1), c(2, NA, 0), c(0, 0, 1), c(1, 0, 0))
  )

  # loci with a different ALT allele can be kept apart (bcftools merge -m
  # none), which makes duplicated positions
  expect_error(
    gt_merge_vcfs(c(vcf_1, vcf_2),
      variants = "union", alt_mismatch = "split",
      backingfile = tempfile(), quiet = TRUE
    ),
    "duplicated loci"
  )
  # (both loci keep their ID, so names are duplicated as well)
  expect_warning(
    expect_warning(
      merged <- gt_merge_vcfs(c(vcf_1, vcf_2),
        variants = "union", alt_mismatch = "split", allow_duplicates = TRUE,
        backingfile = tempfile(), quiet = TRUE
      ),
      "duplicated locus names"
    ),
    "duplicated loci"
  )
  expect_identical(
    show_loci(merged)$position,
    c(100L, 100L, 200L, 300L, 400L, 50L)
  )
  expect_identical(show_loci(merged)$allele_alt[1:2], c("C", "G"))
  expect_equal(
    show_genotypes(merged)[, 1:2],
    cbind(c(1, 2, NA), c(NA, NA, 2))
  )
})

test_that("gt_merge_vcfs matches rbind() and handles REF mismatches", {
  vcf_a <- system.file("extdata", "pop_a.vcf", package = "tidypopgen")
  vcf_b <- system.file("extdata", "pop_b.vcf", package = "tidypopgen")
  # the REF alleles of some loci differ between these files: an error by
  # default (as in bcftools merge)
  expect_error(
    gt_merge_vcfs(c(vcf_a, vcf_b), backingfile = tempfile(), quiet = TRUE),
    "the REF alleles of 7 loci differ between files .* 1:798959 \\(G/A in"
  )
  # which can be dropped
  merged <- gt_merge_vcfs(c(vcf_a, vcf_b),
    ref_mismatch = "drop",
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(
    show_loci(merged)$name,
    c(
      "rs3094315", "rs3131972", "rs12124819", "rs1110052", "rs6657048",
      "rs1240719"
    )
  )
  # or, if the alleles are swapped, merged swapping the genotypes, as rbind()
  # does
  merged <- gt_merge_vcfs(c(vcf_a, vcf_b),
    ref_mismatch = "swap",
    backingfile = tempfile(), quiet = TRUE
  )
  pop_a <- gen_tibble(vcf_a, backingfile = tempfile(), quiet = TRUE)
  pop_b <- gen_tibble(vcf_b, backingfile = tempfile(), quiet = TRUE)
  bound <- rbind(pop_a, pop_b,
    use_position = TRUE,
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(merged$id, bound$id)
  loci_cols <- c("name", "position", "allele_ref", "allele_alt")
  expect_identical(
    show_loci(merged)[, loci_cols],
    show_loci(bound)[, loci_cols]
  )
  # (rbind() keeps the chromosomes of the first gen_tibble, even those whose
  # loci were all dropped; gt_merge_vcfs() only those with merged loci)
  expect_identical(
    as.character(show_loci(merged)$chromosome),
    as.character(show_loci(bound)$chromosome)
  )
  expect_identical(levels(show_loci(merged)$chromosome), "1")
  expect_identical(show_genotypes(merged), show_genotypes(bound))
})

test_that("gt_merge_vcfs handles samples, ploidy, integer names and gz", {
  # the same sample name in two files is an error, unless force_samples
  # renames it as bcftools does
  vcf_3 <- write_merge_vcf(
    samples = c("a", "d"),
    loci = data.frame(
      chrom = c("1", "1"), pos = c(200, 300), id = c("s2", "s3"),
      ref = c("G", "C"), alt = c("T", "G")
    ),
    genotypes = rbind(c("1", "."), c("0", "1"))
  )
  expect_error(
    gt_merge_vcfs(c(vcf_1, vcf_3), backingfile = tempfile(), quiet = TRUE),
    "force_samples = TRUE"
  )
  # (vcf_3 is haploid: its missing genotypes are those of the merged FBM)
  merged <- gt_merge_vcfs(c(vcf_1, vcf_3),
    force_samples = TRUE,
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(merged$id, c("a", "b", "2:a", "d"))
  expect_identical(indiv_ploidy(merged), c(2L, 2L, 1L, 1L))
  expect_equal(show_genotypes(merged), cbind(c(0, 1, 1, NA), c(2, NA, 0, 1)))

  # integer names, a gzipped file and two cores
  gz_file <- tempfile(fileext = ".vcf.gz")
  gz_con <- gzfile(gz_file, "w")
  writeLines(readLines(vcf_2), gz_con)
  close(gz_con)
  merged <- gt_merge_vcfs(c(vcf_1, gz_file),
    variants = "union", names_as_int = TRUE, n_cores = 2,
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(show_loci(merged)$name, 1:4)
  expect_equal(
    show_genotypes(merged),
    cbind(c(0, 1, 1), c(2, NA, 0), c(NA, NA, 1), c(1, 0, NA))
  )

  # files with no loci in common
  vcf_4 <- write_merge_vcf(
    samples = "e",
    loci = data.frame(chrom = "3", pos = 1, id = "x", ref = "A", alt = "C"),
    genotypes = rbind("0/1")
  )
  expect_error(
    gt_merge_vcfs(c(vcf_1, vcf_4), backingfile = tempfile(), quiet = TRUE),
    "no loci in common"
  )
  expect_error(gt_merge_vcfs(vcf_1), "two or more VCF files")
})
