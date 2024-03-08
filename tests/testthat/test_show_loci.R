testthat::test_that("show_loci gets and sets information",{
  test_ind_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,2),
                          c(2,1,0,NA,0,NA),
                          c(2,2,0,0,1,NA))
  test_loci <- data.frame(name=c(paste0("rs",1:4),"x1","x2"),
                          chromosome=c(1,1,1,1,2,2),
                          position=c(3,5,65,343,23,456),
                          allele_ref = c("a","t","c","g","c","t"),
                          allele_alt = c("t","c", NA,"c","g","a"))
  test_gen <- gen_tibble(test_ind_meta, test_genotypes, test_loci)
  # check that we retrieve the info we put in (as a tibble)
  expect_identical(show_loci(test_gen),as_tibble(test_loci))
  # now change it directly on the genotype column
  test_loci2 <- test_loci %>% dplyr::mutate(chromosome = "new")
  show_loci(test_gen$genotypes) <- test_loci2
  expect_identical(show_loci(test_gen), as_tibble(test_loci2))
  test_loci3 <- test_loci %>% dplyr::mutate(chromosome = "newer")
  show_loci(test_gen) <- test_loci3
  expect_identical(show_loci(test_gen), as_tibble(test_loci3))
  # with some proper dplyr
  show_loci(test_gen) <- show_loci(test_gen) %>% mutate(chromosome="old")
  expect_true(all(show_loci(test_gen)$chromosome=="old"))
})
