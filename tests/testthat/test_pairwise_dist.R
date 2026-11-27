# 3 individuals and 6 loci, with missing genotypes
test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, NA, 0, 0),
  c(2, NA, 0, 0, 1, 1)
)
test_indiv_meta <- data.frame(
  id = c("a", "b", "c"),
  population = c("pop1", "pop1", "pop2")
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

# the upper triangle of a matrix (the pairs of a tidy tibble, in its order)
upper <- function(m) unclass(m)[upper.tri(m)]

test_that("pairwise_dist scores pairs of genotypes as ngsDist", {
  d <- pairwise_dist(test_gt, as_matrix = TRUE, correction = "none")
  expect_s3_class(d, "pairwise_matrix")
  expect_identical(dimnames(d), list(test_gt$id, test_gt$id))
  expect_identical(unname(diag(d)), c(0, 0, 0))
  # by hand: a and b share loci 1, 2, 3, 5, 6 (scores 0.5, 0, 0, 0.5, 0), a
  # and c loci 1, 3, 4, 5, 6 (0.5, 0, 0.5, 0, 0.5), b and c loci 1, 3, 5, 6
  # (0, 0, 0.5, 0.5)
  expect_equal(d["a", "b"], 1 / 5)
  expect_equal(d["a", "c"], 1.5 / 5)
  expect_equal(d["b", "c"], 1 / 4)
  expect_identical(attr(d, "n_loci")[upper.tri(d)], c(5, 5, 4))
  # 1 minus the proportion of alleles identical by state
  ibs <- pairwise_ibs(test_gt, as_matrix = TRUE)
  expect_equal(upper(d), 1 - upper(ibs))
  # with the nucleotide score, two heterozygotes score 0.5 (a and b at locus
  # 2, a and c at locus 5): 1 minus the allele sharing
  nuc <- pairwise_dist(test_gt,
    as_matrix = TRUE, score = "nucleotide",
    correction = "none"
  )
  expect_equal(nuc["a", "b"], 1.5 / 5)
  expect_equal(nuc["a", "c"], 2 / 5)
  expect_equal(nuc["b", "c"], 1 / 4)
  expect_equal(
    upper(nuc),
    1 - upper(pairwise_allele_sharing(test_gt, as_matrix = TRUE))
  )
  expect_identical(attr(nuc, "n_loci"), attr(d, "n_loci"))
  expect_identical(unname(diag(nuc)), c(0, 0, 0))
})

test_that("pairwise_dist corrects the distances", {
  d <- pairwise_dist(test_gt, as_matrix = TRUE, correction = "none")
  # the log correction is the default
  expect_identical(
    pairwise_dist(test_gt, as_matrix = TRUE),
    pairwise_dist(test_gt, as_matrix = TRUE, correction = "log")
  )
  expect_equal(
    upper(pairwise_dist(test_gt, as_matrix = TRUE)),
    -log(1 - upper(d))
  )
  expect_equal(
    upper(pairwise_dist(test_gt, as_matrix = TRUE, correction = "jc")),
    -3 / 4 * log(1 - 4 / 3 * upper(d))
  )
  # opposite homozygotes at every locus: an infinite log distance, and no
  # Jukes-Cantor distance; individuals without loci in common: NA
  far_gt <- gen_tibble(
    x = rbind(c(0, 0, NA), c(2, 2, NA), c(NA, NA, 1)),
    loci = test_loci[1:3, ],
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  expect_warning(
    far_log <- pairwise_dist(far_gt, as_matrix = TRUE),
    "infinite"
  )
  expect_identical(upper(far_log), c(Inf, NA, NA))
  expect_warning(
    far_jc <- pairwise_dist(far_gt, as_matrix = TRUE, correction = "jc"),
    "Jukes-Cantor"
  )
  expect_identical(upper(far_jc), c(NA_real_, NA, NA))
  expect_identical(attr(far_jc, "n_loci")[upper.tri(far_jc)], c(2, 0, 0))
})

test_that("pairwise_dist gives a tidy tibble or a matrix", {
  d <- pairwise_dist(test_gt, as_matrix = TRUE)
  d_tidy <- pairwise_dist(test_gt)
  expect_s3_class(d_tidy, "pairwise_tbl")
  expect_identical(names(d_tidy), c("item1", "item2", "value", "n_loci"))
  expect_identical(d_tidy$item1, c("a", "a", "b"))
  expect_identical(d_tidy$item2, c("b", "c", "c"))
  expect_identical(d_tidy$value, upper(d))
  expect_identical(d_tidy$n_loci, attr(d, "n_loci")[upper.tri(d)])
  expect_error(pairwise_dist(test_gt[1, ]), "at least 2 individuals")
  expect_error(pairwise_dist(test_gt, score = "foo"), "should be one of")
})

test_that("pairwise_dist is a generic with methods for gen_tibbles", {
  d <- pairwise_dist(test_gt, as_matrix = TRUE)
  d_tidy <- pairwise_dist(test_gt)
  # the method for the genotypes column labels the individuals with the names
  # of the vector, and the method for gen_tibbles with their ids
  expect_identical(pairwise_dist(test_gt$genotypes, as_matrix = TRUE), d)
  expect_identical(pairwise_dist(test_gt$genotypes), d_tidy)
  renamed_gt <- test_gt %>% mutate(id = toupper(id))
  expect_identical(
    dimnames(pairwise_dist(renamed_gt, as_matrix = TRUE)),
    list(renamed_gt$id, renamed_gt$id)
  )
  expect_identical(pairwise_dist(renamed_gt)$item1, c("A", "A", "B"))
  expect_error(pairwise_dist(test_gt, foo = 1), "must be empty")
  # grouping is ignored, and tibbles that lost their class work
  expect_identical(pairwise_dist(test_gt %>% group_by(population)), d_tidy)
  expect_identical(pairwise_dist(as_tibble(test_gt)), d_tidy)
  # other packages can add methods for their own classes of genotypes, with
  # their own arguments; the numbers of loci are given if the matrix has them
  registerS3method(
    "genotypes_tbl_class", "dummy_genotypes",
    function(genotypes, ...) "dummy_tbl",
    envir = asNamespace("tidypopgen")
  )
  registerS3method(
    "pairwise_dist", "dummy_genotypes",
    function(x, as_matrix = FALSE, prior = "flat", ...) {
      structure(matrix(0.1, length(x), length(x)),
        n_loci = matrix(7, length(x), length(x)),
        class = c("pairwise_matrix", "matrix", "array")
      )
    },
    envir = asNamespace("tidypopgen")
  )
  dummy_tbl <- tibble(
    id = c("a", "b", "c"),
    genotypes = structure(1:3, class = "dummy_genotypes")
  )
  dummy_dist <- pairwise_dist(dummy_tbl, prior = "other")
  expect_s3_class(dummy_dist, "pairwise_tbl")
  expect_identical(dummy_dist$item2, c("b", "c", "c"))
  expect_identical(dummy_dist$value, c(0.1, 0.1, 0.1))
  expect_identical(dummy_dist$n_loci, c(7, 7, 7))
})
