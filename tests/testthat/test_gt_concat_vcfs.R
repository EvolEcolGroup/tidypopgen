# the same three samples, with chromosome 1 in one file and chromosome 2 in
# another
concat_loci_1 <- data.frame(
  chrom = "1", pos = c(100, 200), id = c("s1", "."),
  ref = c("A", "G"), alt = c("C", "T")
)
concat_geno_1 <- rbind(c("0/1", "1/1", "0/0"), c("./.", "0/1", "1/1"))
concat_loci_2 <- data.frame(
  chrom = "2", pos = c(50, 60, 70), id = c("s3", "s4", "s5"),
  ref = c("T", "C", "A"), alt = c("A", ".", "G")
)
concat_geno_2 <- rbind(
  c("0/0", "0/1", "1/1"), c("0/0", "0/0", "0/0"), c("1/1", "1/1", "0/1")
)
concat_vcf_1 <- write_merge_vcf(c("a", "b", "c"), concat_loci_1, concat_geno_1)
concat_vcf_2 <- write_merge_vcf(c("a", "b", "c"), concat_loci_2, concat_geno_2)

test_that("gt_concat_vcfs reads files as gen_tibble reads them concatenated", {
  concat_gt <- gt_concat_vcfs(c(concat_vcf_1, concat_vcf_2),
    names_as_int = FALSE, backingfile = tempfile(), quiet = TRUE
  )
  # the same data in a single VCF
  single_vcf <- write_merge_vcf(c("a", "b", "c"),
    rbind(concat_loci_1, concat_loci_2),
    rbind(concat_geno_1, concat_geno_2)
  )
  single_gt <- gen_tibble(single_vcf, backingfile = tempfile(), quiet = TRUE)
  expect_identical(concat_gt$id, single_gt$id)
  expect_identical(show_loci(concat_gt), show_loci(single_gt))
  expect_identical(show_genotypes(concat_gt), show_genotypes(single_gt))
  expect_identical(show_loci(concat_gt)$name[1:3], c("s1", "1_200", "s3"))
  # the gen_tibble is saved, and can be loaded back
  reloaded <- gt_load(
    bigstatsr::sub_bk(gt_get_file_names(concat_gt)[2], ".gt")
  )
  expect_identical(show_genotypes(reloaded), show_genotypes(concat_gt))

  # by default (unlike gen_tibble()), loci are named with integers
  int_gt <- gt_concat_vcfs(c(concat_vcf_1, concat_vcf_2),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(show_loci(int_gt)$name, 1:5)
  expect_identical(show_genotypes(int_gt), show_genotypes(single_gt))
})

test_that("gt_concat_vcfs reorders the samples of each file", {
  # file 2 has the samples in another order (c, a, b), and is gzipped
  reordered_vcf <- write_merge_vcf(c("c", "a", "b"),
    concat_loci_2, concat_geno_2[, c(3, 1, 2)],
    file = tempfile(fileext = ".vcf")
  )
  gz_file <- tempfile(fileext = ".vcf.gz")
  gz_con <- gzfile(gz_file, "w")
  writeLines(readLines(reordered_vcf), gz_con)
  close(gz_con)
  reordered_gt <- gt_concat_vcfs(c(concat_vcf_1, gz_file),
    n_cores = 2, backingfile = tempfile(), quiet = TRUE
  )
  ordered_gt <- gt_concat_vcfs(c(concat_vcf_1, concat_vcf_2),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(reordered_gt$id, c("a", "b", "c"))
  expect_identical(show_genotypes(reordered_gt), show_genotypes(ordered_gt))
})

test_that("gt_concat_vcfs checks the samples and the loci", {
  # a sample is missing, or there is another sample
  other_vcf <- write_merge_vcf(c("a", "b", "d"), concat_loci_2, concat_geno_2)
  expect_error(
    gt_concat_vcfs(c(concat_vcf_1, other_vcf),
      backingfile = tempfile(), quiet = TRUE
    ),
    paste(
      "should have the same samples, but .* lacks samples of the first file",
      "\\(c\\) and has samples that are not in the first file \\(d\\)"
    )
  )
  # duplicated samples within a file
  dup_vcf <- write_merge_vcf(c("a", "b", "b"), concat_loci_2, concat_geno_2)
  expect_error(
    gt_concat_vcfs(c(concat_vcf_1, dup_vcf),
      backingfile = tempfile(), quiet = TRUE
    ),
    "file 2 has duplicated sample names"
  )
  expect_error(gt_concat_vcfs(concat_vcf_1), "two or more VCF files")
  # the same locus in two files
  expect_error(
    gt_concat_vcfs(
      c(
        concat_vcf_1,
        write_merge_vcf(c("a", "b", "c"), concat_loci_1, concat_geno_1)
      ),
      backingfile = tempfile(), quiet = TRUE
    ),
    "duplicated loci"
  )
  # files with no biallelic SNPs
  indel_loci <- data.frame(
    chrom = "1", pos = 100, id = "i1", ref = "AT", alt = "A"
  )
  indels <- rbind(c("0/1", "0/0", "1/1"))
  expect_error(
    gt_concat_vcfs(
      c(
        write_merge_vcf(c("a", "b", "c"), indel_loci, indels),
        write_merge_vcf(c("a", "b", "c"), indel_loci, indels)
      ),
      backingfile = tempfile(), quiet = TRUE
    ),
    "no biallelic SNPs"
  )
})

test_that("gt_concat_vcfs keeps the largest ploidy of each individual", {
  # a file of the X chromosome, in which b is haploid
  x_vcf <- write_merge_vcf(c("a", "b", "c"),
    data.frame(chrom = "X", pos = 10, id = "x1", ref = "G", alt = "C"),
    rbind(c("0/1", "1", "1/1"))
  )
  msgs <- capture_messages(
    x_gt <- gt_concat_vcfs(c(x_vcf, concat_vcf_1), backingfile = tempfile())
  )
  expect_match(msgs, "concatenating 2 VCF files with 3 samples", all = FALSE)
  expect_identical(indiv_ploidy(x_gt), c(2L, 2L, 2L))
  expect_equal(show_genotypes(x_gt)[, 1], c(1, 1, 2))
  # a haploid individual (a) whose genotypes are all missing in another file,
  # where its ploidy can only be guessed (as that of the other individuals):
  # it stays haploid, whatever the order of the files
  haploid_vcf <- write_merge_vcf(c("a", "b", "c"),
    data.frame(
      chrom = "1", pos = c(10, 20), id = ".", ref = "G", alt = "C"
    ),
    rbind(c("1", "0/1", "1/1"), c("0", "0/0", "0/1"))
  )
  missing_vcf <- write_merge_vcf(c("a", "b", "c"),
    data.frame(chrom = "2", pos = 10, id = ".", ref = "G", alt = "C"),
    rbind(c(".", "0/1", "1/1"))
  )
  file_orders <- list(
    c(haploid_vcf, missing_vcf),
    c(missing_vcf, haploid_vcf)
  )
  for (files in file_orders) {
    concat_gt <- gt_concat_vcfs(files, backingfile = tempfile(), quiet = TRUE)
    expect_identical(indiv_ploidy(concat_gt), c(1L, 2L, 2L))
    chr_1 <- show_loci(concat_gt)$chromosome == "1"
    expect_equal(show_genotypes(concat_gt)[1, chr_1], c(1, 0))
    expect_true(is.na(show_genotypes(concat_gt)[1, !chr_1]))
  }
})
