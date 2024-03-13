testthat::test_that("as_genlight creates a valid object",{
  library(adegenet)
  test_ind_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,0),
                          c(2,1,0,0,0,0),
                          c(2,2,0,0,1,1))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=c(1,1,1,1,2,2),
                          position=c(3,5,65,343,23,456),
                          allele_ref = c("a","t","c","g","c","t"),
                          allele_alt = c("t","c", NA,"c","g","a"))
  test_gen <- gen_tibble(test_ind_meta, test_genotypes, test_loci)
  test_gl <- as_genlight(test_gen)
  expect_true(validObject(test_gl))
})
