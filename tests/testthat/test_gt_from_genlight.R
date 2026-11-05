skip_if_not_installed("adegenet")

test_that("gt can convert from genlight", {
  x <- new("genlight",
    list(
      indiv1 = c(1, 1, 0, 1, 1, 0),
      indiv2 = c(2, 1, 1, 0, 0, 0)
    ),
    ploidy = c(2, 2),
    loc.names = paste0("locus", 1:6),
    chromosome = c("chr1", "chr1", "chr2", "chr2", "chr3", "chr3"),
    position = c(100, 200, 150, 250, 300, 400),
    loc.all = c("A/T", "C/G", "G/C", "A/T", "T/C", "G/A"),
    pop = c("pop1", "pop2")
  )

  file <- paste0(tempfile(), "gt_from_genlight")

  new_gt <- gt_from_genlight(x, backingfile = file)
  expect_true(inherits(new_gt, "gen_tbl"))

  expect_true(all(file.exists(gt_get_file_names(new_gt))))
  # check that gt_get_file_names(new_gt)[1] ends with .rds
  expect_true(grepl(".rds$", gt_get_file_names(new_gt)[1]))
  # and gt_get_file_names(new_gt)[2] ends with .bk
  expect_true(grepl(".bk$", gt_get_file_names(new_gt)[2]))
})

test_that("gt can convert from genlight without pop", {
  x <- new("genlight",
    list(
      indiv1 = c(1, 1, 0, 1, 1, 0),
      indiv2 = c(2, 1, 1, 0, 0, 0)
    ),
    ploidy = c(2, 2),
    loc.names = paste0("locus", 1:6),
    chromosome = c("chr1", "chr1", "chr2", "chr2", "chr3", "chr3"),
    position = c(100, 200, 150, 250, 300, 400),
    loc.all = c("A/T", "C/G", "G/C", "A/T", "T/C", "G/A")
  )
  new_gt <- gt_from_genlight(x, backingfile = tempfile())
  expect_true(inherits(new_gt, "gen_tbl"))
  expect_equal(names(new_gt), c("id", "genotypes"))
  expect_equal(new_gt$id, c("indiv1", "indiv2"))
  expect_equal(
    show_genotypes(new_gt),
    rbind(c(1, 1, 0, 1, 1, 0), c(2, 1, 1, 0, 0, 0))
  )
})

test_that("gt_as_genlight and gt_from_genlight round trip", {
  test_loci <- data.frame(
    name = paste0("rs", 1:4),
    chromosome = paste0("chr", c(1, 1, 2, 2)),
    position = as.integer(c(3, 5, 23, 456)),
    genetic_dist = as.double(rep(0, 4)),
    allele_ref = c("A", "T", "C", "G"),
    allele_alt = c("T", "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = rbind(c(1, 1, 0, 2), c(2, 1, 0, NA), c(2, 2, 0, 0)),
    loci = test_loci,
    indiv_meta = data.frame(id = c("a", "b", "c"), population = "pop1"),
    backingfile = tempfile(),
    quiet = TRUE
  )
  test_gl <- gt_as_genlight(test_gt)
  # genlight counts the second allele, i.e. the alternate allele
  expect_equal(adegenet::alleles(test_gl), c("A/T", "T/C", "C/G", "G/A"))
  expect_equal(
    as.character(adegenet::chromosome(test_gl)),
    test_loci$chromosome
  )
  expect_equal(adegenet::position(test_gl), test_loci$position)
  new_gt <- gt_from_genlight(test_gl, backingfile = tempfile())
  expect_equal(show_genotypes(new_gt), show_genotypes(test_gt))
  expect_equal(show_loci(new_gt), show_loci(test_gt))
  expect_equal(new_gt$id, test_gt$id)
  expect_equal(as.character(new_gt$population), test_gt$population)
  # genlight can not store missing alleles
  show_loci(test_gt)$allele_alt[3] <- NA
  expect_warning(
    test_gl <- gt_as_genlight(test_gt),
    "the alleles \\(loc.all\\) were not stored"
  )
  expect_null(adegenet::alleles(test_gl))
  expect_equal(adegenet::position(test_gl), test_loci$position)
})

test_that("error with non-diploid genlight", {
  x <- new("genlight",
    list(
      indiv1 = c(1, 1, 0, 1, 1, 0),
      indiv2 = c(2, 1, 1, 0, 0, 0)
    ),
    ploidy = c(2, 3),
    loc.names = paste0("locus", 1:6),
    chromosome = c("chr1", "chr1", "chr2", "chr2", "chr3", "chr3"),
    position = c(100, 200, 150, 250, 300, 400),
    loc.all = c("A/T", "C/G", "G/C", "A/T", "T/C", "G/A"),
    pop = c("pop1", "pop2")
  )
  expect_error(
    gt_from_genlight(x),
    "Currently only diploid genlight objects are supported"
  )
})

test_that("error with null slots", {
  x <- new(
    "genlight",
    list(
      indiv1 = c(1, 1, 0, 1, 1, 0),
      indiv2 = c(2, 1, 1, 0, 0, 0)
    )
  )
  expect_error(
    gt_from_genlight(x),
    "The genlight object has one or more required slots that are NULL"
  )
})
