skip_if_not_installed("hierfstat")

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


test_that("pairwise_grm compute correctly", {
  # expect error if the tibble is not grouped
  coa_hier <- hierfstat::beta.dosage(test_genotypes, inb = FALSE)
  grm_hier <- 2 * coa_hier
  grm_gt <- pairwise_grm(test_gt)
  # remove the special class to allow comparison with the reference
  class(grm_gt) <- setdiff(class(grm_gt), "pairwise_matrix")
  expect_true(all.equal(grm_hier, grm_gt, check.attributes = FALSE))
})

test_that("pairwise_grm(method = 'plink') matches plink --make-rel", {
  # 60 individuals with 0 to 90% missing genotypes and a monomorphic locus,
  # and the output of plink --make-rel square, made by the script
  # create_plink_pca_fixture.R in data-raw
  fixture <- test_path("testdata", "plink_pca", "missing")
  missing_gt <- gen_tibble(paste0(fixture, ".bed"),
    backingfile = tempfile("grm_"), quiet = TRUE
  )
  plink_grm <- pairwise_grm(missing_gt, method = "plink", n_cores = 2)
  plink_rel <- unname(as.matrix(read.table(paste0(fixture, ".rel"))))
  # PLINK writes 6 significant digits
  expect_lt(max(abs(unclass(plink_grm) - plink_rel)), 1e-5)
  expect_s3_class(plink_grm, "pairwise_matrix")
  expect_equal(rownames(plink_grm), missing_gt$id)
  # the default estimator is unchanged, and differs from PLINK's
  expect_false(isTRUE(all.equal(
    unclass(pairwise_grm(missing_gt)), unclass(plink_grm),
    check.attributes = FALSE
  )))
  expect_error(
    pairwise_grm(missing_gt, allele_sharing_mat = diag(60), method = "plink"),
    "only used with method"
  )
})
