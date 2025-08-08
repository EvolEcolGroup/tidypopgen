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

  test_genlight <- gt_as_genlight(test_gt)
  expect_true(all.equal(
    show_genotypes(test_gt),
    as.matrix(test_genlight),
    check.attributes = FALSE
  ))
})
