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
