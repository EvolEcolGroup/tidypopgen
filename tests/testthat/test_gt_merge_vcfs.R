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
  # intersection (the default): sites in both files (keeping the VCF IDs as
  # loci names)
  merged <- gt_merge_vcfs(c(vcf_1, vcf_2),
    names_as_int = FALSE,
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
    variants = "union", names_as_int = FALSE,
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
        names_as_int = FALSE, backingfile = tempfile(), quiet = TRUE
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
    ref_mismatch = "drop", names_as_int = FALSE,
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
    ref_mismatch = "swap", names_as_int = FALSE,
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

  # integer names (the default, unlike in gen_tibble()), a gzipped file and
  # two cores
  gz_file <- tempfile(fileext = ".vcf.gz")
  gz_con <- gzfile(gz_file, "w")
  writeLines(readLines(vcf_2), gz_con)
  close(gz_con)
  merged <- gt_merge_vcfs(c(vcf_1, gz_file),
    variants = "union", n_cores = 2,
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

test_that("gt_merge_vcfs checks its arguments", {
  expect_error(
    gt_merge_vcfs(c(vcf_1, tempfile(fileext = ".vcf"))),
    "these VCF files do not exist"
  )
  txt_file <- tempfile(fileext = ".txt")
  file.copy(vcf_1, txt_file)
  expect_error(gt_merge_vcfs(c(vcf_1, txt_file)), "should be VCF files")
  expect_error(gt_merge_vcfs(c(vcf_1, vcf_1)), "more than once")
  expect_error(
    gt_merge_vcfs(c(vcf_1, vcf_2), valid_alleles = c("A", "0")),
    "'0' cannot be a valid allele"
  )
  expect_error(gt_merge_vcfs(c(vcf_1, vcf_2), n_cores = 0), "positive integer")

  # missing_to_ref is only used with the union of the sites; by default, the
  # backing files are next to the first file
  msgs <- capture_messages(
    merged <- gt_merge_vcfs(c(vcf_1, vcf_2), missing_to_ref = TRUE)
  )
  expect_match(msgs, "only used with variants = 'union'", all = FALSE)
  expect_identical(
    basename(gt_get_file_names(merged)[2]),
    paste0(sub("\\.vcf$", "", basename(vcf_1)), "_merged.bk")
  )
})

test_that("gt_merge_vcfs checks the samples and genotypes of each file", {
  loci <- data.frame(
    chrom = "1", pos = c(100, 200), id = c("s1", "s2"),
    ref = c("A", "G"), alt = c("C", "T")
  )
  vcf_a <- write_merge_vcf("a", loci, rbind("0/1", "1/1"))
  two_samples <- rbind(c("0/1", "0/0"), c("1/1", "0/1"))
  # duplicated sample names within a file
  vcf_dup <- write_merge_vcf(c("b", "b"), loci, two_samples)
  expect_error(
    gt_merge_vcfs(c(vcf_a, vcf_dup), backingfile = tempfile(), quiet = TRUE),
    "file 2 has duplicated sample names"
  )
  # renaming duplicated samples ("a" in file 2 becomes "2:a") can not create
  # names that are already used
  vcf_2a <- write_merge_vcf(c("a", "2:a"), loci, two_samples)
  expect_error(
    gt_merge_vcfs(c(vcf_2a, vcf_a),
      force_samples = TRUE, backingfile = tempfile(), quiet = TRUE
    ),
    "created more duplicated names: 2:a"
  )
  # a VCF without samples (sites only)
  vcf_sites <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t100\ts1\tA\tC\t.\tPASS\t."
  ), vcf_sites)
  expect_error(
    gt_merge_vcfs(c(vcf_a, vcf_sites), backingfile = tempfile(), quiet = TRUE),
    "has no samples"
  )
  # a locus without alleles, but with genotypes
  no_alleles <- loci
  no_alleles$ref[1] <- "."
  no_alleles$alt[1] <- "."
  expect_error(
    gt_merge_vcfs(
      c(
        write_merge_vcf("a", no_alleles, rbind("0/1", "1/1")),
        write_merge_vcf("b", no_alleles, rbind("0/0", "0/1"))
      ),
      backingfile = tempfile(), quiet = TRUE
    ),
    "missing both reference and alternate alleles"
  )
  # a variant line with too few genotypes (the backing files are removed)
  vcf_short <- write_merge_vcf(c("b", "c"), loci, two_samples)
  vcf_lines <- readLines(vcf_short)
  n_lines <- length(vcf_lines)
  vcf_lines[n_lines] <- sub("\t[^\t]*$", "", vcf_lines[n_lines])
  writeLines(vcf_lines, vcf_short)
  bk <- tempfile()
  expect_error(
    gt_merge_vcfs(c(vcf_a, vcf_short), backingfile = bk, quiet = TRUE),
    "has fewer genotypes than individuals"
  )
  expect_false(any(file.exists(paste0(bk, c(".bk", ".rds")))))
  # a genotype with more alleles than the ploidy (from the first variant)
  vcf_haploid <- write_merge_vcf("b", loci, rbind("1", "0/1"))
  expect_error(
    gt_merge_vcfs(c(vcf_a, vcf_haploid),
      backingfile = tempfile(), quiet = TRUE
    ),
    "has more than max_ploidy alleles"
  )
})

test_that("gt_merge_vcfs skips loci that are not biallelic SNPs", {
  vcf_a <- write_merge_vcf("a",
    data.frame(
      chrom = "1", pos = c(100, 200), id = c("s1", "s2"),
      ref = c("A", "G"), alt = c("C", "T")
    ),
    genotypes = rbind("0/1", "1/1")
  )
  # an indel between the two SNPs
  vcf_indel <- write_merge_vcf("b",
    data.frame(
      chrom = "1", pos = c(100, 150, 200), id = c("s1", "i1", "s2"),
      ref = c("A", "AT", "G"), alt = c("C", "A", "T")
    ),
    genotypes = rbind("0/1", "1/1", "0/0")
  )
  merged <- gt_merge_vcfs(c(vcf_a, vcf_indel),
    names_as_int = FALSE, backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(show_loci(merged)$name, c("s1", "s2"))
  expect_equal(show_genotypes(merged), rbind(c(1, 2), c(1, 0)))
})

test_that("gt_merge_vcfs merges files with many loci", {
  # (enough positions to grow the hash table of the sites)
  vcf_a <- write_merge_vcf("a",
    data.frame(chrom = "1", pos = 1:800, id = ".", ref = "A", alt = "C"),
    genotypes = matrix("0/1", 800, 1)
  )
  vcf_b <- write_merge_vcf("b",
    data.frame(chrom = "1", pos = 401:1200, id = ".", ref = "A", alt = "C"),
    genotypes = matrix("1/1", 800, 1)
  )
  merged <- gt_merge_vcfs(c(vcf_a, vcf_b),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(show_loci(merged)$position, 401:800)
  merged <- gt_merge_vcfs(c(vcf_a, vcf_b),
    variants = "union", backingfile = tempfile(), quiet = TRUE
  )
  expect_identical(show_loci(merged)$position, 1:1200)
})

test_that("gt_merge_vcfs keeps loci at the same position in a file apart", {
  # a site with three ALT alleles, split over three lines (as with bcftools
  # norm -m-), and a file with one of them
  vcf_split <- write_merge_vcf("a",
    data.frame(
      chrom = "1", pos = c(100, 100, 100, 200),
      id = c("m1", "m2", "m3", "s2"), ref = c("A", "A", "A", "G"),
      alt = c("C", "G", "T", "T")
    ),
    genotypes = rbind("0/1", "0/0", "0/0", "1/1")
  )
  vcf_b <- write_merge_vcf("b",
    data.frame(
      chrom = "1", pos = c(100, 200), id = c("x", "s2"),
      ref = c("A", "G"), alt = c("G", "T")
    ),
    genotypes = rbind("1/1", "0/1")
  )
  merged <- gt_merge_vcfs(c(vcf_split, vcf_b),
    names_as_int = FALSE, backingfile = tempfile(), quiet = TRUE
  )
  # the locus of file b is merged with the line with the same ALT allele
  expect_identical(show_loci(merged)$name, c("m2", "s2"))
  expect_identical(show_loci(merged)$allele_alt, c("G", "T"))
  expect_equal(show_genotypes(merged), rbind(c(0, 2), c(2, 1)))
})

test_that("the C++ merge helpers check their inputs", {
  # (these inputs are always consistent when they come from gt_merge_vcfs()
  # or gt_concat_vcfs())
  expect_error(
    merge_sites_new(split_alt = FALSE, ref_rule = 3L, missing_alleles = "."),
    "ref_rule should be 0, 1 or 2"
  )
  sites <- merge_sites_new(
    split_alt = FALSE, ref_rule = 0L, missing_alleles = "."
  )
  expect_error(
    merge_sites_add_file(sites,
      chromosome = 1L, position = c(1L, 2L),
      allele_ref = "A", allele_alt = "C"
    ),
    "should have a chromosome, a position and two alleles"
  )
  expect_identical(
    merge_sites_add_file(sites,
      chromosome = c(1L, 1L), position = c(1L, 2L),
      allele_ref = c("A", "G"), allele_alt = c("C", ".")
    ),
    1:2
  )
  expect_error(merge_sites_file_cols(sites, 2L, c(0L, 1L)), "not available")
  expect_error(merge_sites_file_cols(sites, 1L, 0L), "one value per site")
  # the sites of a file are released once their columns are used
  expect_identical(
    merge_sites_file_cols(sites, 1L, c(1L, -1L)),
    list(col = c(1L, -1L), swap = logical(0))
  )
  expect_error(merge_sites_file_cols(sites, 1L, c(0L, 1L)), "not available")
  expect_true(is.na(merge_sites_result(sites)$allele_alt[2]))

  # a VCF with one sample and two SNPs, and an FBM with 2 rows and 2 columns
  vcf_path <- write_merge_vcf("a",
    data.frame(
      chrom = "1", pos = c(100, 200), id = c("s1", "s2"),
      ref = c("A", "G"), alt = c("C", "T")
    ),
    genotypes = rbind("0/1", "1/1")
  )
  code256 <- rep(NA_real_, 256)
  code256[1:3] <- 0:2
  fbm <- bigstatsr::FBM.code256(2, 2, code = code256, backingfile = tempfile())
  merge_genotypes <- function(path = vcf_path, row_offset = 0L,
                              row = integer(0), col = c(0L, 1L),
                              swap = logical(0), missing_value = 3L) {
    vcf_merge_genotypes(path,
      BM = fbm, row_offset = row_offset, row = row, col = col, swap = swap,
      ploidy = 2L, missing_value = missing_value
    )
  }
  merge_genotypes(row = 1L)
  expect_equal(fbm[, ], rbind(c(0, 0), c(1, 2)))
  expect_error(merge_genotypes(row_offset = 2L), "does not have enough rows")
  expect_error(merge_genotypes(row = c(0L, 1L)), "one value per individual")
  expect_error(merge_genotypes(row = 2L), "'row' is out of the bounds")
  expect_error(merge_genotypes(col = c(0L, 2L)), "'col' is out of the bounds")
  expect_error(merge_genotypes(swap = TRUE), "one value per biallelic snp")
  expect_error(merge_genotypes(missing_value = 0L), "between 1 and 255")
  # the file must have as many SNPs as when its loci were read
  expect_error(merge_genotypes(col = 0L), "has more biallelic SNPs than")
  expect_error(merge_genotypes(col = c(0L, 1L, -1L)), "has fewer biallelic")
  short_path <- tempfile(fileext = ".vcf")
  writeLines(c(readLines(vcf_path), "1\t300\ts3\tA"), short_path)
  expect_error(merge_genotypes(path = short_path), "fewer than 5 fields")
})
