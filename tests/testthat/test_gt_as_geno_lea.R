test_that("gt_as_geno_lea errors on polyploid data", {
  test_indiv_meta <- data.frame(id = letters[1:4], population = rep("popA", 4))
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
  tetra_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    ploidy = 4,
    quiet = TRUE
  )
  expect_error(
    gt_as_geno_lea(tetra_gt, file = paste0(tempfile(), ".geno")),
    "this function only works on diploid data"
  )
})

test_that("gt_as_geno_lea writes a .geno file for diploid data", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, NA),
    c(2, 1, 0, 0),
    c(2, 2, 0, 1)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:4),
    chromosome = c(1, 1, 1, 2),
    position = c(3, 5, 65, 456),
    genetic_dist = as.double(rep(0, 4)),
    allele_ref = c("A", "T", "C", "G"),
    allele_alt = c("T", "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  geno_file <- gt_as_geno_lea(test_gt, file = paste0(tempfile(), ".geno"))
  # the .geno format has one row per locus, with 1 character per individual
  # concatenated with no separator (9 for missing)
  geno_lines <- readLines(geno_file)
  geno_matrix <- do.call(rbind, strsplit(geno_lines, ""))
  expected <- t(test_genotypes)
  expected[is.na(expected)] <- 9
  storage.mode(expected) <- "character"
  dimnames(geno_matrix) <- NULL
  dimnames(expected) <- NULL
  expect_equal(geno_matrix, expected)
})
