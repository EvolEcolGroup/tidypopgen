# a small gen_tibble with autosomes, sex chromosomes and mtDNA, named in
# several ways
test_indiv_meta <- data.frame(
  id = c("a", "b", "c"),
  population = c("pop1", "pop1", "pop2")
)
test_genotypes <- rbind(
  c(1, 1, 0, 2, 2, 0, 2, 0),
  c(2, 1, 0, 1, 0, NA, 0, 2),
  c(2, 2, 0, 0, 2, 2, 2, 0)
)
test_loci <- data.frame(
  name = paste0("rs", 1:8),
  chromosome = c("1", "1", "23", "chrX", "chrX", "Y", "chrM", "PAR1"),
  position = as.integer(c(3, 5, 65, 1000, 3e6, 23, 456, 10)),
  genetic_dist = as.double(rep(0, 8)),
  allele_ref = c("A", "T", "C", "G", "C", "T", "A", "G"),
  allele_alt = c("T", "C", "G", "C", "G", "A", "G", "A")
)

test_that("chromosome types are recognised from names, not numbers", {
  test_gt <- gen_tibble(
    x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta,
    quiet = TRUE, backingfile = tempfile()
  )
  expect_identical(
    as.character(loci_chromosome_types(test_gt)),
    c("autosome", "autosome", "autosome", "X", "X", "Y", "MT", "PAR")
  )
  expect_identical(
    levels(loci_chromosome_types(test_gt$genotypes)),
    c("autosome", "X", "Y", "Z", "W", "PAR", "MT")
  )
  # (sorted in the C locale, as the order of levels depends on the locale)
  types <- gt_chromosome_types(test_gt) %>%
    dplyr::arrange(chromosome, type)
  expect_identical(types$chromosome, c("1", "23", "PAR1", "Y", "chrM", "chrX"))
  expect_identical(
    as.character(types$type),
    c("autosome", "autosome", "PAR", "Y", "MT", "X")
  )
  expect_identical(types$n_loci, c(2L, 1L, 1L, 1L, 1L, 2L))
  # the types are found in a loci table subset
  test_gt <- test_gt %>%
    select_loci_if(loci_chromosome_types(genotypes) == "autosome")
  expect_identical(loci_names(test_gt), c("rs1", "rs2", "rs3"))
})

test_that("gt_set_chromosome_types declares types and PARs", {
  test_gt <- gen_tibble(
    x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta,
    quiet = TRUE, backingfile = tempfile()
  )
  # the human codes of PLINK
  human_gt <- gt_set_chromosome_types(test_gt, chr_set = "human")
  expect_identical(
    as.character(loci_chromosome_types(human_gt))[3],
    "X"
  )
  # explicit types take precedence, and bring chromosomes back to autosomes
  decl_gt <- gt_set_chromosome_types(test_gt,
    types = c("23" = "MT", chrX = "autosome"), chr_set = 22
  )
  expect_identical(
    as.character(loci_chromosome_types(decl_gt))[3:5],
    c("MT", "autosome", "autosome")
  )
  # a call without a declaration clears it
  expect_null(attr(gt_set_chromosome_types(decl_gt)$genotypes, "chr_types"))
  # pseudo-autosomal regions of a build (rs4 at 1000 is in PAR1 of b37)
  par_gt <- gt_set_chromosome_types(test_gt, par = "b37")
  expect_identical(
    as.character(loci_chromosome_types(par_gt))[4:5],
    c("PAR", "X")
  )
  par_types <- gt_chromosome_types(par_gt) %>%
    dplyr::filter(chromosome == "chrX")
  expect_identical(as.character(par_types$type), c("X", "PAR"))
  expect_identical(par_types$n_loci, c(1L, 1L))
  # or given for a chromosome (positions included)
  par_gt <- gt_set_chromosome_types(test_gt,
    par = data.frame(chromosome = "chrX", start = 1000, end = 3e6)
  )
  expect_identical(
    as.character(loci_chromosome_types(par_gt))[4:5],
    c("PAR", "PAR")
  )
  # bad declarations
  expect_error(
    gt_set_chromosome_types(test_gt, types = c("X", "Y")),
    "named character vector"
  )
  expect_error(
    gt_set_chromosome_types(test_gt, types = c("23" = "sex")),
    "should be one of"
  )
  expect_error(gt_set_chromosome_types(test_gt, chr_set = "dog"), "chr_set")
  expect_error(gt_set_chromosome_types(test_gt, par = "b39"), "'par'")
})

test_that("declared chromosome types are kept with the gen_tibble", {
  test_gt <- gen_tibble(
    x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta,
    quiet = TRUE, backingfile = tempfile()
  ) %>%
    gt_set_chromosome_types(chr_set = "human", par = "b37")
  expected <- attr(test_gt$genotypes, "chr_types")
  # subsetting individuals and loci
  sub_gt <- test_gt %>%
    dplyr::filter(id != "b") %>%
    select_loci(-1)
  expect_identical(attr(sub_gt$genotypes, "chr_types"), expected)
  # saving and loading
  gt_file <- gt_save(test_gt, quiet = TRUE)[1]
  expect_identical(
    attr(gt_load(gt_file)$genotypes, "chr_types"),
    expected
  )
  # new backing files
  new_gt <- gt_update_backingfile(sub_gt,
    backingfile = tempfile(),
    quiet = TRUE
  )
  expect_identical(attr(new_gt$genotypes, "chr_types"), expected)
  # binding rows: the union of the declarations, which must agree
  other_gt <- gen_tibble(
    x = test_genotypes, loci = test_loci,
    indiv_meta = data.frame(id = c("d", "e", "f"), population = "pop3"),
    quiet = TRUE, backingfile = tempfile()
  ) %>%
    gt_set_chromosome_types(types = c(chrM = "MT", "26" = "MT"))
  bound_gt <- rbind(test_gt, other_gt,
    quiet = TRUE,
    backingfile = tempfile()
  )
  expect_identical(
    attr(bound_gt$genotypes, "chr_types")$types,
    c(expected$types, "chrM" = "MT")
  )
  other_gt <- gt_set_chromosome_types(other_gt, types = c("23" = "autosome"))
  expect_error(
    rbind(test_gt, other_gt, quiet = TRUE, backingfile = tempfile()),
    "different types"
  )
})

test_that("chromosome types are found for integer chromosomes", {
  vcf_path <- write_merge_vcf(c("a", "b"),
    data.frame(
      chrom = c("chr1", "chrX", "chrM"), pos = c(10, 20, 30),
      id = c("rs1", "rs2", "rs3"), ref = c("A", "C", "G"),
      alt = c("T", "G", "A")
    ),
    rbind(c("0/1", "1/1"), c("0/1", "1/1"), c("0/0", "1/1"))
  )
  test_gt <- gen_tibble(vcf_path,
    backingfile = tempfile(), quiet = TRUE,
    chromosomes_as_int = TRUE
  )
  expect_true(is.integer(show_loci(test_gt)$chromosome))
  expect_identical(
    as.character(loci_chromosome_types(test_gt)),
    c("autosome", "X", "MT")
  )
  expect_identical(
    gt_chromosome_types(test_gt)$chromosome,
    c("chr1", "chrX", "chrM")
  )
})
