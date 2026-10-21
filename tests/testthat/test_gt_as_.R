skip_if_not_installed("adegenet")
skip_if_not_installed("hierfstat")

test_that("show_loci gets and sets information", {
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, NA, 0, 0),
    c(2, NA, 0, 0, 1, 1),
    c(1, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(0, 0, 0, 0, NA, 1),
    c(0, 1, 1, 0, 1, NA)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e", "f", "g"),
    population = c("pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3")
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(3, 5, 65, 343, 23, 456)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  test_hier <- gt_as_hierfstat(test_gt)
  test_genind <- gt_as_genind(test_gt)
  # now test that the two objects are identical
  # (minus the row names, which are abritrary)
  expect_true(all.equal(
    test_hier,
    hierfstat::genind2hierfstat(test_genind),
    check.attributes = FALSE
  ))

  # genlight can not store the missing allele of rs3
  expect_warning(
    test_genlight <- gt_as_genlight(test_gt),
    "the alleles \\(loc.all\\) were not stored"
  )
  expect_true(all.equal(
    show_genotypes(test_gt),
    as.matrix(test_genlight),
    check.attributes = FALSE
  ))
})

# a tetraploid gen_tibble, used to check that hierfstat (which does not
# support polyploid data) errors cleanly, and that genind (which does)
# decodes dosage into allele counts correctly
tetraploid_gt <- function() {
  test_indiv_meta <- data.frame(
    id = letters[1:4],
    population = c("popA", "popA", "popB", "popB")
  )
  test_genotypes <- rbind(
    c(1, 4, 0, 2),
    c(4, 4, 0, 0),
    c(0, 2, 0, 4),
    c(2, 0, 4, NA)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:4),
    chromosome = c(1, 1, 2, 2),
    position = c(3, 5, 23, 456),
    genetic_dist = as.double(rep(0, 4)),
    allele_ref = c("A", "T", "C", "G"),
    allele_alt = c("T", "C", "G", "A")
  )
  gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    ploidy = 4,
    quiet = TRUE
  )
}

# a,b diploid; c,d tetraploid
mixed_ploidy_gt <- function() {
  test_indiv_meta <- data.frame(id = letters[1:4], population = rep("popA", 4))
  test_genotypes <- rbind(
    c(1, 2, 0, 1),
    c(0, 2, 1, 0),
    c(4, 4, 0, 2),
    c(0, 1, 4, 4)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:4),
    chromosome = c(1, 1, 2, 2),
    position = c(3, 5, 23, 456),
    genetic_dist = as.double(rep(0, 4)),
    allele_ref = c("A", "T", "C", "G"),
    allele_alt = c("T", "C", "G", "A")
  )
  gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    ploidy = c(2, 2, 4, 4),
    quiet = TRUE
  )
}

test_that("gt_as_hierfstat errors on polyploid data", {
  expect_error(
    gt_as_hierfstat(tetraploid_gt()),
    "this function only works on diploid or pseudohaploid data"
  )
  expect_error(
    gt_as_hierfstat(mixed_ploidy_gt()),
    "this function only works on diploid or pseudohaploid data"
  )
})

test_that("gt_as_hierfstat keeps population labels and loci names", {
  test_gt <- gen_tibble(
    x = rbind(c(0, 1, 2), c(1, 2, NA), c(2, 0, 0)),
    loci = data.frame(
      name = c("rs1", "rs2", "rs3"),
      chromosome = 1,
      position = 1:3,
      genetic_dist = 0,
      allele_ref = c("A", "C", "G"),
      allele_alt = c("G", "T", "A")
    ),
    # population labels that look like genotype codes
    indiv_meta = data.frame(
      id = c("a", "b", "c"),
      population = c("1", "2", "2")
    ),
    quiet = TRUE
  )
  test_hier <- gt_as_hierfstat(test_gt)
  expect_equal(test_hier$pop, factor(c("1", "2", "2")))
  expect_equal(names(test_hier), c("pop", "rs1", "rs2", "rs3"))
  expect_equal(test_hier$rs1, c(11, 12, 22))
  expect_equal(test_hier$rs2, c(12, 22, 11))
  expect_equal(test_hier$rs3, c(22, NA, 11))
})

test_that("gt_as_hierfstat treats pseudohaploid dosage as diploid", {
  test_genotypes <- rbind(
    c(2, NA, 0, NA, 2, 0),
    c(2, 2, 0, NA, 0, NA),
    c(2, NA, 0, 0, 1, 1),
    c(2, 0, 0, 1, 0, 0)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d"),
    population = rep("pop1", 4)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(3, 5, 65, 343, 23, 456)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt <- gt_pseudohaploid(test_gt)
  expect_true(is_pseudohaploid(test_gt))
  test_hier <- gt_as_hierfstat(test_gt)
  # pseudohaploid individuals a and b are homozygous diploids
  expect_equal(unname(unlist(test_hier[1, -1])), c(22, NA, 11, NA, 22, 11))
  expect_equal(unname(unlist(test_hier[2, -1])), c(22, 22, 11, NA, 11, NA))
})

test_that("gt_as_genind works on uniform polyploid data", {
  test_gt <- tetraploid_gt()
  test_genind <- gt_as_genind(test_gt)
  expect_equal(unname(adegenet::ploidy(test_genind)), rep(4L, 4))

  # individual a, locus rs1: dosage 1 out of ploidy 4
  allele_tab <- adegenet::tab(test_genind)
  expect_equal(unname(allele_tab["a", "rs1.2"]), 1)
  expect_equal(unname(allele_tab["a", "rs1.1"]), 3)
  # individual d, locus rs4: missing dosage stays missing
  expect_true(all(is.na(allele_tab["d", c("rs4.1", "rs4.2")])))
})

test_that("gt_as_genind works on mixed-ploidy data", {
  test_gt <- mixed_ploidy_gt()
  test_genind <- gt_as_genind(test_gt)
  expect_equal(unname(adegenet::ploidy(test_genind)), c(2L, 2L, 4L, 4L))

  allele_tab <- adegenet::tab(test_genind)
  # individual c (tetraploid), locus rs4: dosage 2 out of ploidy 4
  expect_equal(unname(allele_tab["c", "rs4.1"]), 2)
  expect_equal(unname(allele_tab["c", "rs4.2"]), 2)
  # individual a (diploid), locus rs1: dosage 1 out of ploidy 2
  expect_equal(unname(allele_tab["a", "rs1.1"]), 1)
  expect_equal(unname(allele_tab["a", "rs1.2"]), 1)
})

test_that("gt_as_genind replaces '.' in loci names", {
  test_gt <- gen_tibble(
    x = rbind(c(0, 1, 2), c(1, NA, 1)),
    loci = data.frame(
      name = c("NC_031965.1_10", "NC_031965.1_20", "chr2_5"),
      chromosome = "NC_031965.1",
      position = 1:3,
      genetic_dist = 0,
      allele_ref = c("A", "C", "G"),
      allele_alt = c("G", "T", "T")
    ),
    indiv_meta = data.frame(id = c("a", "b"), population = "pop1"),
    quiet = TRUE
  )
  # genind uses '.' to separate loci and alleles
  expect_warning(
    test_genind <- gt_as_genind(test_gt),
    "character '.' detected in names of loci; replacing with '_'"
  )
  expect_equal(
    unname(adegenet::locNames(test_genind)),
    c("NC_031965_1_10", "NC_031965_1_20", "chr2_5")
  )
  allele_tab <- adegenet::tab(test_genind)
  expect_equal(unname(allele_tab["a", "NC_031965_1_20.2"]), 1)
  expect_equal(unname(allele_tab["b", "chr2_5.2"]), 1)
})

test_that("gt_as_genind treats pseudohaploid dosage as diploid", {
  test_genotypes <- rbind(
    c(2, NA, 0, NA, 2, 0),
    c(2, 2, 0, NA, 0, NA),
    c(2, NA, 0, 0, 1, 1),
    c(2, 0, 0, 1, 0, 0)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d"),
    population = rep("pop1", 4)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(3, 5, 65, 343, 23, 456)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt <- gt_pseudohaploid(test_gt)
  expect_equal(indiv_ploidy(test_gt), c(1, 1, 2, 2))

  test_genind <- gt_as_genind(test_gt)
  # pseudohaploid individuals are represented as ploidy 2, matching their
  # doubled dosage storage (0/2, never 1)
  expect_equal(unname(adegenet::ploidy(test_genind)), c(2L, 2L, 2L, 2L))
  allele_tab <- adegenet::tab(test_genind)
  expect_equal(unname(allele_tab["a", "rs1.2"]), 2)
})
