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
  test_fbm <- tidypopgen:::.gt_get_fbm(test_gt)
  test_king <- snp_king(test_fbm)
  # king by hand
  test_king_r <- king_r(show_genotypes(test_gt))
  expect_identical(test_king_r, test_king)
  # check that we get the same result if we split the operation into two blocks
  test_king_2blocks <- snp_king(test_fbm, block.size = 3)
  expect_identical(test_king_2blocks, test_king)

  # now estimate it with gen_tibble
  test_king_gt <- pairwise_king(test_gt, as_matrix = TRUE)
  # remove the special class to allow comparison with the reference
  class(test_king_gt) <- setdiff(class(test_king_gt), "pairwise_matrix")
  expect_true(all.equal(test_king, test_king_gt, check.attributes = FALSE))

  # now test with missing data
  test_na_gt <- gen_tibble(
    system.file("extdata/related/families.bed", package = "tidypopgen"),
    quiet = TRUE,
    backingfile = tempfile(),
    valid_alleles = c("1", "2")
  )
  test_na_fbm <- tidypopgen:::.gt_get_fbm(test_na_gt)
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
  families_fbm <- tidypopgen:::.gt_get_fbm(families)
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

test_that("snp_king is correct across bit-packing boundaries", {
  # more than 64 individuals and loci, with missing data, so that genotypes
  # are packed into several words and blocks of individuals
  set.seed(123)
  geno_vals <- sample(c(0:2, NA), 70 * 150,
    replace = TRUE,
    prob = c(0.4, 0.3, 0.25, 0.05)
  )
  geno_mat <- matrix(geno_vals, nrow = 70)
  big_gt <- gen_tibble(
    x = geno_mat,
    loci = data.frame(
      name = paste0("rs", 1:150), chromosome = 1L, position = 1:150,
      genetic_dist = 0, allele_ref = "A", allele_alt = "T"
    ),
    indiv_meta = data.frame(id = paste0("i", 1:70), population = "p1"),
    quiet = TRUE
  )
  test_fbm <- tidypopgen:::.gt_get_fbm(big_gt)
  # reordered subsets of individuals and loci, and multiple cores
  rows <- c(70:3, 1)
  cols <- c(150:80, 2:60)
  test_king <- snp_king(test_fbm,
    ind.row = rows, ind.col = cols,
    block.size = 70, n_cores = 2
  )
  expect_equal(test_king, king_r(geno_mat[rows, cols]), tolerance = 1e-12)
})

test_that("snp_king decodes all bytes of the FBM through its code256", {
  # bytes decoded as genotypes (0 to 2, imputed 4 to 6, dosages 7, 107 and
  # 207), and as missing (3, and the dosage 0.09 of byte 16), in blocks of
  # individuals with and without bytes above 15
  # (the code of bigsnpr's dosages, CODE_DOSAGE)
  code <- c(0, 1, 2, NA, 0, 1, 2, seq(0, 2, by = 0.01), rep(NA, 48))
  set.seed(7)
  byte_vals <- sample(c(0:7, 16, 107, 207), 40 * 150,
    replace = TRUE, prob = c(rep(0.12, 8), 0.02, 0.01, 0.01)
  )
  bytes <- matrix(byte_vals, nrow = 40)
  # the first 16 individuals only have bytes below 16
  bytes[1:16, ][bytes[1:16, ] > 15] <- 1
  test_fbm <- bigstatsr::FBM.code256(40, 150,
    code = code, init = bytes, backingfile = tempfile()
  )
  geno_mat <- test_fbm[]
  geno_mat[!geno_mat %in% 0:2] <- NA
  expect_equal(snp_king(test_fbm, n_cores = 2), king_r(geno_mat),
    tolerance = 1e-12
  )
})

test_that("pairwise_king uses imputed genotypes only if set to", {
  geno_na <- test_genotypes
  geno_na[2, 4] <- NA
  geno_na[3, 2] <- NA
  gt_na <- gen_tibble(
    x = geno_na, loci = test_loci, indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  king_na <- pairwise_king(gt_na, as_matrix = TRUE)
  gt_imputed <- gt_impute_simple(gt_na, method = "mode")
  # imputed genotypes are missing unless the gen_tibble is set to use them
  gt_set_imputed(gt_imputed, set = FALSE)
  expect_identical(pairwise_king(gt_imputed, as_matrix = TRUE), king_na)
  # when used, they are the same as genotypes that were never missing
  gt_set_imputed(gt_imputed, set = TRUE)
  gt_filled <- gen_tibble(
    x = show_genotypes(gt_imputed), loci = test_loci,
    indiv_meta = test_indiv_meta, quiet = TRUE
  )
  expect_identical(
    pairwise_king(gt_imputed, as_matrix = TRUE),
    pairwise_king(gt_filled, as_matrix = TRUE)
  )
})

test_that("pairwise_king is a generic with methods for gen_tibbles", {
  king_mat <- pairwise_king(test_gt, as_matrix = TRUE)
  king_tidy <- pairwise_king(test_gt)
  # the method for the genotypes column labels the individuals with the names
  # of the vector, and the method for gen_tibbles with their ids
  expect_identical(pairwise_king(test_gt$genotypes, as_matrix = TRUE), king_mat)
  expect_identical(pairwise_king(test_gt$genotypes), king_tidy)
  renamed_gt <- test_gt %>% mutate(id = toupper(id))
  expect_identical(
    dimnames(pairwise_king(renamed_gt, as_matrix = TRUE)),
    list(renamed_gt$id, renamed_gt$id)
  )
  expect_identical(pairwise_king(renamed_gt)$item1, c("A", "A", "B"))
  # arguments are matched by the methods, also by position
  expect_identical(pairwise_king(test_gt, TRUE, 2, 1), king_mat)
  expect_error(pairwise_king(test_gt, foo = 1), "must be empty")
  # tibbles that lost their class work
  expect_identical(pairwise_king(as_tibble(test_gt)), king_tidy)
  # other packages can add methods for their own classes of genotypes
  registerS3method(
    "genotypes_tbl_class", "dummy_genotypes",
    function(genotypes, ...) "dummy_tbl",
    envir = asNamespace("tidypopgen")
  )
  registerS3method(
    "pairwise_king", "dummy_genotypes",
    function(x, as_matrix = FALSE, ...) {
      structure(diag(0.5, length(x)),
        class = c("pairwise_matrix", "matrix", "array")
      )
    },
    envir = asNamespace("tidypopgen")
  )
  dummy_tbl <- tibble(
    id = c("a", "b", "c"),
    genotypes = structure(1:3, class = "dummy_genotypes")
  )
  dummy_king <- pairwise_king(dummy_tbl)
  expect_s3_class(dummy_king, "pairwise_tbl")
  expect_identical(dummy_king$item2, c("b", "c", "c"))
  expect_identical(dummy_king$value, c(0, 0, 0))
  # grouped tibbles of other classes of genotypes are computed group by group
  dummy_grouped <- dummy_tbl %>%
    mutate(pop = c("p1", "p2", "p1")) %>%
    group_by(pop)
  dummy_king <- pairwise_king(dummy_grouped)
  expect_identical(dummy_king$group, "p1")
  expect_identical(dummy_king$item2, "c")
  expect_length(pairwise_king(dummy_grouped, as_matrix = TRUE), 2)
})

test_that("pairwise_king on a grouped gen_tibble works within groups", {
  # more than 64 individuals and loci, with missing data, in interleaved
  # groups, so that genotypes are packed into several words and blocks of
  # individuals; one group of a single individual, and an empty group
  set.seed(42)
  geno_vals <- sample(c(0:2, NA), 140 * 150,
    replace = TRUE,
    prob = c(0.4, 0.3, 0.25, 0.05)
  )
  pops <- sample(c("A", "B", "C"), 140, replace = TRUE, prob = c(0.6, 0.2, 0.2))
  pops[17] <- "D"
  big_gt <- gen_tibble(
    x = matrix(geno_vals, nrow = 140),
    loci = data.frame(
      name = paste0("rs", 1:150), chromosome = 1L, position = 1:150,
      genetic_dist = 0, allele_ref = "A", allele_alt = "T"
    ),
    indiv_meta = data.frame(
      id = paste0("i", 1:140),
      population = factor(pops, levels = c("C", "A", "E", "D", "B")),
      sex = rep(c("F", "M"), 70)
    ),
    quiet = TRUE
  )
  # a subset of the loci, in reversed order
  big_gt <- big_gt %>% select_loci(c(150:80, 2:60))
  grouped <- big_gt %>% group_by(population, .drop = FALSE)
  expected <- group_map(grouped, ~ pairwise_king(.x, as_matrix = TRUE))
  # a list of matrices, one per group (empty for the empty group), with one
  # block of loci or several, and several cores
  expect_identical(pairwise_king(grouped, as_matrix = TRUE), expected)
  expect_identical(
    pairwise_king(grouped, as_matrix = TRUE, block_size = 70, n_cores = 2),
    expected
  )
  expect_identical(dim(expected[[3]]), c(0L, 0L))
  # the tidy version has the pairs of all groups, after a group column that
  # keeps the type of the grouping variable
  king_tidy <- pairwise_king(grouped)
  expect_named(king_tidy, c("group", "item1", "item2", "value"))
  expect_identical(levels(king_tidy$group), levels(big_gt$population))
  # (a plain tibble: a pairwise_tbl holds all the pairs of a set of items)
  expect_identical(
    king_tidy[-1],
    tibble::as_tibble(dplyr::bind_rows(lapply(expected, function(king_g) {
      class(king_g) <- NULL
      tidypopgen:::tidy_dist_matrix(king_g)
    })))
  )
  group_sizes <- vapply(expected, nrow, integer(1))
  expect_identical(
    as.character(king_tidy$group),
    rep(levels(big_gt$population), group_sizes * (group_sizes - 1) / 2)
  )
  # several grouping variables give one column each
  grouped_2 <- big_gt %>% group_by(population, sex)
  king_tidy_2 <- pairwise_king(grouped_2)
  expect_named(king_tidy_2, c("population", "sex", "item1", "item2", "value"))
  expect_identical(
    pairwise_king(grouped_2, as_matrix = TRUE),
    group_map(grouped_2, ~ pairwise_king(.x, as_matrix = TRUE))
  )
  # a grouped tibble that has lost the class of its gen_tibble
  expect_identical(
    pairwise_king(as_tibble(big_gt) %>% group_by(population), as_matrix = TRUE),
    group_map(
      big_gt %>% group_by(population),
      ~ pairwise_king(.x, as_matrix = TRUE)
    )
  )
  expect_error(pairwise_king(grouped, foo = 1), "must be empty")
})
