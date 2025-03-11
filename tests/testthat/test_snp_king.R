# create file
test_indiv_meta <- data.frame(
  id = c("a", "b", "c"),
  population = c("pop1", "pop1", "pop2")
)
test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, 0, 0, 0),
  c(2, 2, 0, 0, 1, 1)
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


# function to compute robust king in R
king_r <- function(X_mat) { # nolint start
  X_mat0 <- X_mat == 0
  X_mat0[is.na(X_mat0)] <- 0
  X_mat1 <- X_mat == 1
  X_mat1[is.na(X_mat1)] <- 0
  X_mat2 <- X_mat == 2
  X_mat2[is.na(X_mat2)] <- 0
  king_num <- (X_mat1 %*%
    t(X_mat1) -
    2 * ((X_mat0) %*% t(X_mat2) + (X_mat2) %*% t(X_mat0)))
  X_mat_valid <- !is.na(X_mat)
  N_mat_Aa_i <- X_mat1 %*% t(X_mat_valid)
  N_mat_Aa_j <- t(N_mat_Aa_i)
  king_num /
    (2 * pmin(N_mat_Aa_i, N_mat_Aa_j)) +
    0.5 -
    0.25 * (N_mat_Aa_i + N_mat_Aa_j) / pmin(N_mat_Aa_i, N_mat_Aa_j)
} # nolint end


# this also tests show_genotypes and show_loci
test_that("snp_king and pairwise_king compute king-robust correctly", {
  test_fbm <- tidypopgen:::.gt_get_bigsnp(test_gt)$genotypes
  test_king <- snp_king(test_fbm)
  # king by hand
  test_king_r <- king_r(show_genotypes(test_gt))
  expect_identical(test_king_r, test_king)
  # check that we get the same result if we split the operation into two blocks
  test_king_2blocks <- snp_king(test_fbm, block.size = 3)
  expect_identical(test_king_2blocks, test_king)

  # now estimate it with gen_tibble
  test_king_gt <- pairwise_king(test_gt, as_matrix = TRUE)
  expect_true(all.equal(test_king, test_king_gt, check.attributes = FALSE))

  # now test with missing data
  test_na_gt <- gen_tibble(
    system.file("extdata/related/families.bed", package = "tidypopgen"),
    quiet = TRUE,
    backingfile = tempfile(),
    valid_alleles = c("1", "2")
  )
  test_na_fbm <- tidypopgen:::.gt_get_bigsnp(test_na_gt)$genotypes
  test_na_king <- snp_king(test_na_fbm)
  # king by hand
  test_na_king_r <- king_r(show_genotypes(test_na_gt))
  expect_identical(test_na_king_r, test_na_king)

  # check that we get the same result if we split the operation into two blocks
  test_na_king_2blocks <- snp_king(test_na_fbm, block.size = 300)
  expect_identical(test_na_king_2blocks, test_na_king)
})

test_that("snp_king gives the same results as plink", {
  # Create gentibble for our data
  bed_path <- system.file(
    "extdata/related/families.bed",
    package = "tidypopgen"
  )
  families_bigsnp_path <- bigsnpr::snp_readBed(
    bed_path,
    backingfile = tempfile()
  )
  families <- gen_tibble(
    families_bigsnp_path,
    quiet = TRUE,
    valid_alleles = c("1", "2")
  )

  # Get snp_king results
  families_fbm <- tidypopgen:::.gt_get_bigsnp(families)$genotypes
  families_king <- snp_king(families_fbm)

  # Read in results from king -b families_k.bed --kinship
  king <- read.table(
    system.file("extdata/related/test_king.kin0", package = "tidypopgen"),
    header = FALSE
  )

  # Create empty matrix
  king_matrix <- matrix(nrow = 12, ncol = 12)

  # Fill matrix
  for (i in seq_len(nrow(king))) {
    x <- as.numeric(king[i, 2])
    y <- as.numeric(king[i, 3])

    king_matrix[x, y] <- king[i, 8]
    king_matrix[y, x] <- king[i, 8]
  }

  # Replace diagonal
  diag(king_matrix) <- 0.5000

  diff <- king_matrix - families_king

  result <- all.equal(king_matrix, families_king, tolerance = 0.06)

  # Check results are the same
  expect_true(result)
})
