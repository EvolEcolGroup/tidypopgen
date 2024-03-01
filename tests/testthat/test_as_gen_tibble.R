testthat::test_that("as_gen_tibble works on genlight",{
library(adegenet)
test_genlight <- new("genlight", list(indiv1=c(1,1,0,1,1,0), indiv2=c(2,1,1,0,0,0),
                          toto=c(2,2,0,0,4,4)))
test_gen <- as_gen_tibble(test_genlight)
# check that genotypes are preserved
expect_true (all.equal(as.matrix(test_genlight),
                       show_genotypes(test_gen), check.attributes = FALSE))
# check that individual names (the only thing that was present) were preserved
expect_true (all.equal(test_genlight$ind.names,
                       test_gen$id, check.attributes = FALSE))
# now check the error for empty slots
expect_error (as_gen_tibble(test_genlight,ignore_null_slots = FALSE),
              "pop is null ")
})
