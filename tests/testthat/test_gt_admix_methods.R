# skip if admixture is not installed
skip_if(
  (system2("which", args = "admixture", stdout = NULL) != 0) &&
    !requireNamespace("tidygenclust", quietly = TRUE)
)
# create gentibble
vcf_path <-
  system.file(
    "/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
    package = "tidypopgen"
  )
anole_gt <- gen_tibble(
  vcf_path,
  quiet = TRUE,
  backingfile = tempfile("anolis_")
)
pops_path <- system.file(
  "/extdata/anolis/plot_order_punctatus_n46.csv",
  package = "tidypopgen"
)
pops <- readr::read_csv(pops_path, show_col_types = FALSE)
anole_gt <- anole_gt %>% mutate(id = gsub("punc_", "", .data$id, ))
anole_gt <- anole_gt %>% mutate(population = pops$pop[match(pops$ID, .data$id)])

# group by population
anole_gt <- anole_gt %>% group_by(population)

# create reference gt_admix
admix1 <- gt_admixture(
  anole_gt,
  k = 2:3,
  n_runs = 2,
  crossval = TRUE,
  n_cores = 1,
  seed = c(123, 234),
  conda_env = "none"
)

# test combining identical objects
test_that("combine identical gt_admix objects correctly", {
  # create identical admix object
  admix2 <- admix1
  # combine objects, match_attributes = TRUE
  comb1 <- c(admix1, admix2)
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb1, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb1$k, c(
    admix1$k[1], admix1$k[2], admix2$k[1], admix2$k[2],
    admix1$k[3], admix1$k[4], admix2$k[3], admix2$k[4]
  ))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb1$Q) == length(admix1$Q) + length(admix2$Q))
  # check Q matrices match expectation
  expect_equal(comb1$Q[[5]], admix1$Q[[3]])
  expect_equal(comb1$Q[[7]], admix2$Q[[3]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb1$P) == length(admix1$P) + length(admix2$P))
  # check P matrices match expectation
  expect_equal(comb1$P[[2]], admix1$P[[2]])
  expect_equal(comb1$P[[4]], admix2$P[[2]])
  # check logs match expectation
  expect_equal(comb1$log[[6]], admix1$log[[4]])
  expect_equal(comb1$log[[3]], admix2$log[[1]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb1$id, admix1$id)
  expect_equal(comb1$id, admix2$id)
  # check that the combined object has the same group as the original objects
  expect_equal(comb1$group, admix1$group)
  expect_equal(comb1$group, admix2$group)
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb1$algorithm, admix1$algorithm)
  expect_equal(comb1$algorithm, admix2$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb1$cv, c(
    admix1$cv[1], admix1$cv[2], admix2$cv[1],
    admix2$cv[2], admix1$cv[3], admix1$cv[4],
    admix2$cv[3], admix2$cv[4]
  ))
  # check the loglik values are combined correctly
  expect_equal(comb1$loglik, c(
    admix1$loglik[1], admix1$loglik[2],
    admix2$loglik[1], admix2$loglik[2],
    admix1$loglik[3], admix1$loglik[4],
    admix2$loglik[3], admix2$loglik[4]
  ))

  # combine objects, match_attributes = FALSE
  comb2 <- c(admix1, admix2, match_attributes = FALSE)
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb2, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb2$k, c(
    admix1$k[1], admix1$k[2], admix2$k[1], admix2$k[2],
    admix1$k[3], admix1$k[4], admix2$k[3], admix2$k[4]
  ))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb2$Q) == length(admix1$Q) + length(admix2$Q))
  # check Q matrices match expectation
  expect_equal(comb2$Q[[5]], admix1$Q[[3]])
  expect_equal(comb2$Q[[7]], admix2$Q[[3]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb2$P) == length(admix1$P) + length(admix2$P))
  # check P matrices match expectation
  expect_equal(comb2$P[[2]], admix1$P[[2]])
  expect_equal(comb2$P[[4]], admix2$P[[2]])
  # check logs match expectation
  expect_equal(comb2$log[[5]], admix1$log[[3]])
  expect_equal(comb2$log[[3]], admix2$log[[1]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb2$id, admix1$id)
  expect_equal(comb2$id, admix2$id)
  # check that the combined object has the same group as the original objects
  expect_equal(comb2$group, admix1$group)
  expect_equal(comb2$group, admix2$group)
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb2$algorithm, admix1$algorithm)
  expect_equal(comb2$algorithm, admix2$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb2$cv, c(
    admix1$cv[1], admix1$cv[2], admix2$cv[1],
    admix2$cv[2], admix1$cv[3], admix1$cv[4],
    admix2$cv[3], admix2$cv[4]
  ))
  # check the loglik values are combined correctly
  expect_equal(comb2$loglik, c(
    admix1$loglik[1], admix1$loglik[2],
    admix2$loglik[1], admix2$loglik[2],
    admix1$loglik[3], admix1$loglik[4],
    admix2$loglik[3], admix2$loglik[4]
  ))
})

# test combining objects with same attributes but different k and runs
test_that("combine gt_admix objects with different kand runs", {
  # create admix object with different k range
  admix2 <- gt_admixture(
    anole_gt,
    k = 5:6,
    n_runs = 2,
    crossval = TRUE,
    n_cores = 1,
    seed = c(345, 456),
    conda_env = "none"
  )
  # combine objects, match_attributes = TRUE
  comb1 <- c(admix1, admix2)
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb1, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb1$k, c(admix1$k, admix2$k))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb1$Q) == length(admix1$Q) + length(admix2$Q))
  # check Q matrices match expectation
  expect_equal(comb1$Q[[3]], admix1$Q[[3]])
  expect_equal(comb1$Q[[7]], admix2$Q[[3]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb1$P) == length(admix1$P) + length(admix2$P))
  # check P matrices match expectation
  expect_equal(comb1$P[[2]], admix1$P[[2]])
  expect_equal(comb1$P[[6]], admix2$P[[2]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb1$id, admix1$id)
  expect_equal(comb1$id, admix2$id)
  # check that the combined object has the same group as the original objects
  expect_equal(comb1$group, admix1$group)
  expect_equal(comb1$group, admix2$group)
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb1$algorithm, admix1$algorithm)
  expect_equal(comb1$algorithm, admix2$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb1$cv, c(admix1$cv, admix2$cv))
  # check the loglik values are combined correctly
  expect_equal(comb1$loglik, c(admix1$loglik, admix2$loglik))

  # create admix object with different size of k range
  admix3 <- gt_admixture(
    anole_gt,
    k = 3:5,
    n_runs = 2,
    crossval = TRUE,
    n_cores = 1,
    seed = c(345, 456),
    conda_env = "none"
  )
  # combine objects, match_attributes = TRUE
  comb2 <- c(admix1, admix3)
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb2, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb2$k, c(admix1$k, admix3$k))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb2$Q) == length(admix1$Q) + length(admix3$Q))
  # check Q matrices match expectation
  expect_equal(comb2$Q[[4]], admix1$Q[[4]])
  expect_equal(comb2$Q[[9]], admix3$Q[[5]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb2$P) == length(admix1$P) + length(admix3$P))
  # check P matrices match expectation
  expect_equal(comb2$P[[3]], admix1$P[[3]])
  expect_equal(comb2$P[[8]], admix3$P[[4]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb2$id, admix1$id)
  expect_equal(comb2$id, admix3$id)
  # check that the combined object has the same group as the original objects
  expect_equal(comb2$group, admix1$group)
  expect_equal(comb2$group, admix3$group)
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb2$algorithm, admix1$algorithm)
  expect_equal(comb2$algorithm, admix3$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb2$cv, c(admix1$cv, admix3$cv))
  # check the loglik values are combined correctly
  expect_equal(comb2$loglik, c(admix1$loglik, admix3$loglik))

  # same k different no. of runs
  admix4 <- c(admix1, admix1) # 4 runs
  # combine objects, match_attributes = TRUE
  comb3 <- c(admix1, admix4)
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb3, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb3$k, c(
    admix1$k[1], admix1$k[2], admix4$k[1], admix4$k[2],
    admix4$k[3], admix4$k[4], admix1$k[3], admix1$k[4], admix4$k[5],
    admix4$k[6], admix4$k[7], admix4$k[8]
  ))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb3$Q) == length(admix1$Q) + length(admix4$Q))
  # check Q matrices match expectation
  expect_equal(comb3$Q[[2]], admix1$Q[[2]])
  expect_equal(comb3$Q[[9]], admix4$Q[[5]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb3$P) == length(admix1$P) + length(admix4$P))
  # check P matrices match expectation
  expect_equal(comb3$P[[7]], admix1$P[[3]])
  expect_equal(comb3$P[[6]], admix4$P[[4]])
  # check logs match expectation
  expect_equal(comb3$log[[8]], admix1$log[[4]])
  expect_equal(comb3$log[[11]], admix4$log[[7]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb3$id, admix1$id)
  expect_equal(comb3$id, admix4$id)
  # check that the combined object has the same group as the original objects
  expect_equal(comb3$group, admix1$group)
  expect_equal(comb3$group, admix4$group)
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb3$algorithm, admix1$algorithm)
  expect_equal(comb3$algorithm, admix4$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb3$cv, c(
    admix1$cv[1], admix1$cv[2], admix4$cv[1],
    admix4$cv[2], admix4$cv[3], admix4$cv[4],
    admix1$cv[3], admix1$cv[4], admix4$cv[5],
    admix4$cv[6], admix4$cv[7], admix4$cv[8]
  ))
  # check the loglik values are combined correctly
  expect_equal(comb3$loglik, c(
    admix1$loglik[1], admix1$loglik[2],
    admix4$loglik[1], admix4$loglik[2],
    admix4$loglik[3], admix4$loglik[4],
    admix1$loglik[3], admix1$loglik[4],
    admix4$loglik[5], admix4$loglik[6],
    admix4$loglik[7], admix4$loglik[8]
  ))

  # different k and different no. of runs
  admix5 <- gt_admixture(
    anole_gt,
    k = 3:5,
    n_runs = 1,
    crossval = TRUE,
    n_cores = 1,
    seed = 123,
    conda_env = "none"
  )
  # combine objects, match_attributes = TRUE
  comb4 <- c(admix1, admix5)
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb4, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb4$k, c(admix1$k, admix5$k))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb4$Q) == length(admix1$Q) + length(admix5$Q))
  # check Q matrices match expectation
  expect_equal(comb4$Q[[2]], admix1$Q[[2]])
  expect_equal(comb4$Q[[7]], admix5$Q[[3]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb4$P) == length(admix1$P) + length(admix5$P))
  # check P matrices match expectation
  expect_equal(comb4$P[[3]], admix1$P[[3]])
  expect_equal(comb4$P[[6]], admix5$P[[2]])
  # check logs match expectation
  expect_equal(comb4$log[[3]], admix1$log[[3]])
  expect_equal(comb4$log[[5]], admix5$log[[1]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb4$id, admix1$id)
  expect_equal(comb4$id, admix5$id)
  # check that the combined object has the same group as the original objects
  expect_equal(comb4$group, admix1$group)
  expect_equal(comb4$group, admix5$group)
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb4$algorithm, admix1$algorithm)
  expect_equal(comb4$algorithm, admix5$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb4$cv, c(admix1$cv, admix5$cv))
  # check the loglik values are combined correctly
  expect_equal(comb4$loglik, c(admix1$loglik, admix5$loglik))
})

# test non-matching no. individuals, id order or gt_admix class
test_that("gt_admix objects with different individuals arent combined", {
  # same ids but different order in q-matrices
  admix2 <- admix1
  admix2$id[1] <- "GN71"
  admix2$id[2] <- "BM288"
  # combine objects, match_attributes = TRUE
  expect_error(
    c(admix1, admix2),
    paste(
      "All 'gt_admix' objects have the same individual ids, but individuals",
      "are in a different order. Please re-order individuals before",
      "combining objects using gt_admix_reorder_q()"
    )
  )
  # combine objects, match_attributes = FALSE
  expect_error(
    c(admix1, admix2, match_attributes = FALSE),
    paste(
      "All 'gt_admix' objects have the same individual ids, but individuals",
      "are in a different order. Please re-order individuals before",
      "combining objects using gt_admix_reorder_q()"
    )
  )

  # different number of individuals in q-matrices
  admix3 <- admix1
  admix3$Q <- lapply(admix3$Q, function(mat) mat[1:45, , drop = FALSE])
  admix3$id <- admix3$id[1:45]
  admix3$group <- admix3$group[1:45]
  # combine objects, match_attributes = TRUE
  expect_error(
    c(admix1, admix3),
    "All Q matrices must have the same number of individuals \\(rows\\)"
  )
  # combine objects, match_attributes = FALSE
  expect_error(
    c(admix1, admix3, match_attributes = FALSE),
    "All Q matrices must have the same number of individuals \\(rows\\)"
  )

  # non gt_admix class
  admix4 <- admix1
  class(admix4) <- "list"
  # combine objects, match_attributes = TRUE
  expect_error(c(admix1, admix4), "All the objects must be of class gt_admix")
  # combine objects, match_attributes = FALSE
  expect_error(
    c(admix1, admix4, match_attributes = FALSE),
    "All the objects must be of class gt_admix"
  )
})

# test combining objects with different attributes
test_that("combine gt_admix objects with different attributes", {
  admix2 <- gt_admixture(
    anole_gt,
    k = 3:4,
    n_runs = 2,
    crossval = TRUE,
    n_cores = 1,
    seed = c(123, 234),
    conda_env = "none"
  )
  # create admix object with different id
  admix2$id <- paste0(admix2$id, "_diff")
  # combine objects, match_attributes = TRUE
  expect_error(c(admix1, admix2), "Id information does not match")
  # combine objects, match_attributes = FALSE
  expect_warning(
    comb1 <- c(admix1, admix2, match_attributes = FALSE),
    "Id information is present but does not match. Omitting id"
  )
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb1, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb1$k, c(admix1$k, admix2$k))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb1$Q) == length(admix1$Q) + length(admix2$Q))
  # check Q matrices match expectation
  expect_equal(comb1$Q[[3]], admix1$Q[[3]])
  expect_equal(comb1$Q[[7]], admix2$Q[[3]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb1$P) == length(admix1$P) + length(admix2$P))
  # check P matrices match expectation
  expect_equal(comb1$P[[2]], admix1$P[[2]])
  expect_equal(comb1$P[[5]], admix2$P[[1]])
  # check logs match expectation
  expect_equal(comb1$log[[2]], admix1$log[[2]])
  expect_equal(comb1$log[[5]], admix2$log[[1]])
  # check that the combined object has no id
  expect_true(!"id" %in% names(comb1))
  # check that the combined object has the same group as the original objects
  expect_equal(comb1$group, admix1$group)
  expect_equal(comb1$group, admix2$group)
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb1$algorithm, admix1$algorithm)
  expect_equal(comb1$algorithm, admix2$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb1$cv, c(admix1$cv, admix2$cv))
  # check the loglik values are combined correctly
  expect_equal(comb1$loglik, c(admix1$loglik, admix2$loglik))

  # create admix object with different groups
  admix3 <- admix2
  admix3$id <- admix1$id
  admix3$group <- paste0(admix3$group, "_diff")
  # combine objects, match_attributes = TRUE
  expect_error(c(admix1, admix3), "Group information does not match")
  # combine objects, match_attributes = FALSE
  expect_warning(
    comb2 <- c(admix1, admix3, match_attributes = FALSE),
    "Group information is present but does not match. Omitting group"
  )
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb2, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb2$k, c(admix1$k, admix3$k))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb2$Q) == length(admix1$Q) + length(admix3$Q))
  # check Q matrices match expectation
  expect_equal(comb2$Q[[3]], admix1$Q[[3]])
  expect_equal(comb2$Q[[7]], admix3$Q[[3]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb2$P) == length(admix1$P) + length(admix2$P))
  # check P matrices match expectation
  expect_equal(comb2$P[[2]], admix1$P[[2]])
  expect_equal(comb2$P[[6]], admix3$P[[2]])
  # check logs match expectation
  expect_equal(comb2$log[[1]], admix1$log[[1]])
  expect_equal(comb2$log[[8]], admix3$log[[4]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb2$id, admix1$id)
  expect_equal(comb2$id, admix3$id)
  # check that the combined object has no group
  expect_true(!"group" %in% names(comb2))
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb2$algorithm, admix1$algorithm)
  expect_equal(comb2$algorithm, admix3$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb2$cv, c(admix1$cv, admix3$cv))
  # check the loglik values are combined correctly
  expect_equal(comb2$loglik, c(admix1$loglik, admix3$loglik))

  # create admix object with different algorithm
  admix4 <- admix2
  admix4$id <- admix1$id
  admix4$algorithm <- "fastmixture"
  # combine objects, match_attributes = TRUE
  expect_error(c(admix1, admix4), "Algorithm information does not match")
  # combine objects, match_attributes = FALSE
  expect_warning(
    comb3 <- c(admix1, admix4, match_attributes = FALSE),
    "Algorithm information is present but does not match. Omitting algorithm"
  )
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb3, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb3$k, c(admix1$k, admix4$k))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb3$Q) == length(admix1$Q) + length(admix4$Q))
  # check Q matrices match expectation
  expect_equal(comb3$Q[[2]], admix1$Q[[2]])
  expect_equal(comb3$Q[[7]], admix4$Q[[3]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb3$P) == length(admix1$P) + length(admix4$P))
  # check P matrices match expectation
  expect_equal(comb3$P[[3]], admix1$P[[3]])
  expect_equal(comb3$P[[8]], admix4$P[[4]])
  # check logs match expectation
  expect_equal(comb3$log[[3]], admix1$log[[3]])
  expect_equal(comb3$log[[7]], admix4$log[[3]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb3$id, admix1$id)
  expect_equal(comb3$id, admix4$id)
  # check that the combined object has the same group as the original objects
  expect_equal(comb3$group, admix1$group)
  expect_equal(comb3$group, admix4$group)
  # check that the combined object has no algorithm
  expect_true(!"algorithm" %in% names(comb3))
  # check the cv values are combined correctly
  expect_equal(comb3$cv, c(admix1$cv, admix4$cv))
  # check the loglik values are combined correctly
  expect_equal(comb3$loglik, c(admix1$loglik, admix4$loglik))
})


# test combining objects where one has attributes and one doesnt
test_that("combine gt_admix where one had attribute and one doesnt", {
  # one gt_admix has id one doesn't
  admix2 <- admix1
  admix3 <- admix1
  admix3$id <- NULL
  # combine objects, match_attributes = TRUE
  expect_error(c(admix1, admix3), "Id information does not match")
  # combine objects, match_attributes = FALSE
  expect_warning(
    c(admix1, admix3, match_attributes = FALSE),
    "Only some gt_admix objects have id information; using this id"
  )
  # check same no matter order
  expect_error(c(admix2, admix3), "Id information does not match")
  expect_warning(
    c(admix2, admix3, match_attributes = FALSE),
    "Only some gt_admix objects have id information; using this id"
  )
  # check same if combine 3 objects
  expect_error(c(admix1, admix3, admix2), "Id information does not match")
  expect_warning(
    comb1 <- c(admix1, admix3, admix2, match_attributes = FALSE),
    "Only some gt_admix objects have id information; using this id"
  )
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb1, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb1$k, c(
    admix1$k[1], admix1$k[2], admix3$k[1], admix3$k[2],
    admix2$k[1], admix2$k[2], admix1$k[3], admix1$k[4],
    admix3$k[3], admix3$k[4], admix2$k[3], admix2$k[4]
  ))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb1$Q) == length(admix1$Q) + length(admix2$Q)
    + length(admix3$Q)) # nolint
  # check Q matrices match expectation
  expect_equal(comb1$Q[[2]], admix1$Q[[2]])
  expect_equal(comb1$Q[[6]], admix2$Q[[2]])
  expect_equal(comb1$Q[[3]], admix3$Q[[1]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb1$P) == length(admix1$P) + length(admix2$P)
    + length(admix3$P)) # nolint
  # check P matrices match expectation
  expect_equal(comb1$P[[7]], admix1$P[[3]])
  expect_equal(comb1$P[[12]], admix2$P[[4]])
  expect_equal(comb1$P[[4]], admix3$P[[2]])
  # check logs match expectation
  expect_equal(comb1$log[[8]], admix1$log[[4]])
  expect_equal(comb1$log[[12]], admix2$log[[4]])
  expect_equal(comb1$log[[9]], admix3$log[[3]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb1$id, admix1$id)
  expect_equal(comb1$id, admix2$id)
  expect_false(identical(comb1$id, admix3$id))
  # check that the combined object has the same group as the original objects
  expect_equal(comb1$group, admix1$group)
  expect_equal(comb1$group, admix2$group)
  expect_equal(comb1$group, admix3$group)
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb1$algorithm, admix1$algorithm)
  expect_equal(comb1$algorithm, admix2$algorithm)
  expect_equal(comb1$algorithm, admix3$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb1$cv, c(
    admix1$cv[1], admix1$cv[2], admix3$cv[1],
    admix3$cv[2], admix2$cv[1], admix2$cv[2],
    admix1$cv[3], admix1$cv[4], admix3$cv[3],
    admix3$cv[4], admix2$cv[3], admix2$cv[4]
  ))
  # check the loglik values are combined correctly
  expect_equal(comb1$loglik, c(
    admix1$loglik[1], admix1$loglik[2],
    admix3$loglik[1], admix3$loglik[2],
    admix2$loglik[1], admix2$loglik[2],
    admix1$loglik[3], admix1$loglik[4],
    admix3$loglik[3], admix3$loglik[4],
    admix2$loglik[3], admix2$loglik[4]
  ))

  # one gt_admix has group one doesn't
  admix4 <- admix1
  admix4$group <- NULL
  # combine objects, match_attributes = TRUE
  expect_error(c(admix1, admix4), "Group information does not match")
  # combine objects, match_attributes = FALSE
  expect_warning(
    c(admix1, admix4, match_attributes = FALSE),
    "Only some gt_admix objects have group information; using this group"
  )
  # check same no matter order
  expect_error(c(admix4, admix1), "Group information does not match")
  expect_warning(
    c(admix4, admix1, match_attributes = FALSE),
    "Only some gt_admix objects have group information; using this group"
  )
  # check same if combine 3 objects
  expect_error(c(admix1, admix4, admix2), "Group information does not match")
  expect_warning(
    comb2 <- c(admix1, admix2, admix4, match_attributes = FALSE),
    "Only some gt_admix objects have group information; using this group"
  )
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb2, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb1$k, c(
    admix1$k[1], admix1$k[2], admix2$k[1], admix2$k[2],
    admix4$k[1], admix4$k[2], admix1$k[3], admix1$k[4],
    admix2$k[3], admix2$k[4], admix4$k[3], admix4$k[4]
  ))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb2$Q) == length(admix1$Q) + length(admix4$Q) +
    length(admix2$Q)) # nolint
  # check Q matrices match expectation
  expect_equal(comb2$Q[[7]], admix1$Q[[3]])
  expect_equal(comb2$Q[[11]], admix4$Q[[3]])
  expect_equal(comb2$Q[[3]], admix2$Q[[1]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb2$P) == length(admix1$P) + length(admix4$P) +
    length(admix2$P)) # nolint
  # check P matrices match expectation
  expect_equal(comb2$P[[1]], admix1$P[[1]])
  expect_equal(comb2$P[[5]], admix4$P[[1]])
  expect_equal(comb2$P[[8]], admix2$P[[4]])
  # check logs match expectation
  expect_equal(comb2$log[[2]], admix1$log[[2]])
  expect_equal(comb2$log[[6]], admix4$log[[2]])
  expect_equal(comb2$log[[7]], admix2$log[[3]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb2$id, admix1$id)
  expect_equal(comb2$id, admix4$id)
  expect_equal(comb2$id, admix2$id)
  # check that the combined object has correct group
  expect_equal(comb2$group, admix1$group)
  expect_false(identical(comb2$group, admix4$group))
  expect_equal(comb2$group, admix2$group)
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb2$algorithm, admix1$algorithm)
  expect_equal(comb2$algorithm, admix2$algorithm)
  expect_equal(comb2$algorithm, admix4$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb2$cv, c(
    admix1$cv[1], admix1$cv[2], admix2$cv[1],
    admix2$cv[2], admix4$cv[1], admix4$cv[2],
    admix1$cv[3], admix1$cv[4], admix2$cv[3],
    admix2$cv[4], admix4$cv[3], admix4$cv[4]
  ))
  # check the loglik values are combined correctly
  expect_equal(comb2$loglik, c(
    admix1$loglik[1], admix1$loglik[2],
    admix2$loglik[1], admix2$loglik[2],
    admix4$loglik[1], admix4$loglik[2],
    admix1$loglik[3], admix1$loglik[4],
    admix2$loglik[3], admix2$loglik[4],
    admix4$loglik[3], admix4$loglik[4]
  ))

  # one gt_admix has algorithm one doesn't
  admix5 <- admix1
  admix5$algorithm <- NULL
  # combine objects, match_attributes = TRUE
  expect_error(c(admix1, admix5), "Algorithm information does not match")
  # combine objects, match_attributes = FALSE
  expect_warning(
    c(admix1, admix5, match_attributes = FALSE),
    paste(
      "Only some gt_admix objects have algorithm information;",
      "using this algorithm"
    )
  )
  # check same no matter order
  expect_error(c(admix5, admix1), "Algorithm information does not match")
  expect_warning(
    c(admix5, admix1, match_attributes = FALSE),
    paste(
      "Only some gt_admix objects have algorithm information;",
      "using this algorithm"
    )
  )
  # check same if combine 3 objects
  expect_error(
    c(admix1, admix5, admix2),
    "Algorithm information does not match"
  )
  expect_warning(
    comb3 <- c(admix5, admix2, admix1, match_attributes = FALSE),
    paste(
      "Only some gt_admix objects have algorithm information;",
      "using this algorithm"
    )
  )
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb3, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb1$k, c(
    admix5$k[1], admix5$k[2], admix2$k[1], admix2$k[2],
    admix1$k[1], admix1$k[2], admix5$k[3], admix5$k[4],
    admix2$k[3], admix2$k[4], admix2$k[3], admix2$k[4]
  ))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb3$Q) == length(admix1$Q) + length(admix5$Q) +
    length(admix2$Q)) # nolint
  # check Q matrices match expectation
  expect_equal(comb3$Q[[6]], admix1$Q[[2]])
  expect_equal(comb3$Q[[1]], admix5$Q[[1]])
  expect_equal(comb3$Q[[10]], admix2$Q[[4]])
  # check that the combined object has the correct number of P matrices
  expect_true(length(comb3$P) == length(admix1$P) + length(admix5$P) +
    length(admix2$P)) # nolint
  # check P matrices match expectation
  expect_equal(comb3$P[[11]], admix1$P[[3]])
  expect_equal(comb3$P[[7]], admix5$P[[3]])
  expect_equal(comb3$P[[9]], admix2$P[[3]])
  # check logs match expectation
  expect_equal(comb3$log[[12]], admix1$log[[4]])
  expect_equal(comb3$log[[8]], admix5$log[[4]])
  expect_equal(comb3$log[[10]], admix2$log[[4]])
  # check that the combined object has the same id as the original objects
  expect_equal(comb3$id, admix1$id)
  expect_equal(comb3$id, admix5$id)
  expect_equal(comb3$id, admix2$id)
  # check that the combined object has correct group
  expect_equal(comb3$group, admix1$group)
  expect_equal(comb3$group, admix5$group)
  expect_equal(comb3$group, admix2$group)
  # check that combined object has the same algorithm as the original objects
  expect_equal(comb3$algorithm, admix1$algorithm)
  expect_false(identical(comb3$algorithm, admix5$algorithm))
  expect_equal(comb3$algorithm, admix2$algorithm)
  # check the cv values are combined correctly
  expect_equal(comb3$cv, c(
    admix1$cv[1], admix1$cv[2], admix5$cv[1],
    admix5$cv[2], admix2$cv[1], admix2$cv[2],
    admix1$cv[3], admix1$cv[4], admix5$cv[3],
    admix5$cv[4], admix2$cv[3], admix2$cv[4]
  ))
  # check the loglik values are combined correctly
  expect_equal(comb3$loglik, c(
    admix1$loglik[1], admix1$loglik[2],
    admix5$loglik[1], admix5$loglik[2],
    admix2$loglik[1], admix2$loglik[2],
    admix1$loglik[3], admix1$loglik[4],
    admix5$loglik[3], admix5$loglik[4],
    admix2$loglik[3], admix2$loglik[4]
  ))
})

# test combining objects where one has optional extras (cv, Pmatrices, loglik,
# logs ) one doesnt
test_that("combine gt_admix where one has optional extras and one doesn't", {
  # one gt_admix has cv one doesn't
  admix2 <- admix1
  admix3 <- admix1
  admix3$cv <- NULL
  # combine objects, match_attributes = TRUE
  comb1 <- c(admix1, admix2, admix3)
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb1, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb1$k, c(
    admix1$k[1], admix1$k[2], admix2$k[1], admix2$k[2],
    admix3$k[1], admix3$k[2], admix1$k[3], admix1$k[4],
    admix2$k[3], admix2$k[4], admix3$k[3], admix3$k[4]
  ))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb1$Q) == length(admix1$Q) + length(admix2$Q) +
    length(admix3$Q)) # nolint
  # check the combined object doesnt have cv, but has others
  expect_true(!"cv" %in% names(comb1))
  expect_true(
    all(
      c(
        "k", "Q", "P", "log", "loglik", "id",
        "group", "algorithm"
      ) %in% names(comb1)
    )
  )
  # check the same with match_attributes = FALSE
  comb2 <- c(admix1, admix2, admix3, match_attributes = FALSE)
  expect_equal(comb1, comb2)

  # one gt_admix has P-matrices one doesn't
  admix4 <- admix1
  admix4$P <- NULL
  # combine objects, match_attributes = TRUE
  comb3 <- c(admix1, admix4, admix2)
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb3, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb1$k, c(
    admix1$k[1], admix1$k[2], admix4$k[1], admix4$k[2],
    admix2$k[1], admix2$k[2], admix1$k[3], admix1$k[4],
    admix4$k[3], admix4$k[4], admix2$k[3], admix2$k[4]
  ))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb3$Q) == length(admix1$Q) + length(admix4$Q) +
    length(admix2$Q)) # nolint
  # check the combined object doesnt have P-matrices, but has others
  expect_true(!"P" %in% names(comb3))
  expect_true(
    all(
      c(
        "k", "Q", "cv", "log", "loglik", "id", "group",
        "algorithm"
      ) %in% names(comb3)
    )
  )
  # check the same with match_attributes = FALSE
  comb4 <- c(admix1, admix4, admix2, match_attributes = FALSE)
  expect_equal(comb3, comb4)

  # one gt_admix has logs one doesn't
  admix5 <- admix1
  admix5$log <- NULL
  # combine objects, match_attributes = TRUE
  comb5 <- c(admix5, admix1, admix2)
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb5, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb1$k, c(
    admix5$k[1], admix5$k[2], admix1$k[1], admix1$k[2],
    admix2$k[1], admix2$k[2], admix5$k[3], admix5$k[4],
    admix1$k[3], admix1$k[4], admix2$k[3], admix2$k[4]
  ))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb5$Q) == length(admix5$Q) + length(admix1$Q) +
    length(admix2$Q)) # nolint
  # check the combined object doesnt have logs, but has others
  expect_true(!"^log$" %in% names(comb5))
  expect_true(
    all(
      c(
        "k", "Q", "P", "cv", "loglik", "id", "group",
        "algorithm"
      ) %in% names(comb5)
    )
  )
  # check the same with match_attributes = FALSE
  comb6 <- c(admix5, admix1, admix2, match_attributes = FALSE)
  expect_equal(comb5, comb6)

  # one gt_admix has loglik one doesn't
  admix6 <- admix1
  admix6$loglik <- NULL
  # combine objects, match_attributes = TRUE
  comb7 <- c(admix1, admix2, admix6)
  # check that the combined object is of class gt_admix
  expect_true(inherits(comb7, "gt_admix"))
  # check that the combined object has the same k values as the original objects
  expect_equal(comb1$k, c(
    admix1$k[1], admix1$k[2], admix2$k[1], admix2$k[2],
    admix6$k[1], admix6$k[2], admix1$k[3], admix1$k[4],
    admix2$k[3], admix2$k[4], admix6$k[3], admix6$k[4]
  ))
  # check that the combined object has the correct number of Q matrices
  expect_true(length(comb7$Q) == length(admix1$Q) + length(admix2$Q) +
    length(admix6$Q)) # nolint
  # check the combined object doesnt have loglik, but has others
  expect_true(!"loglik" %in% names(comb7))
  expect_true(
    all(
      c(
        "k", "Q", "P", "cv", "log", "id", "group",
        "algorithm"
      ) %in% names(comb7)
    )
  )
  # check the same with match_attributes = FALSE
  comb8 <- c(admix1, admix2, admix6, match_attributes = FALSE)
  expect_equal(comb7, comb8)
})
