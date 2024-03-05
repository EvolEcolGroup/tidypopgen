# this also tests show_genotypes and show_loci
testthat::test_that("gen_tibble stores data correctly",{
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
  expect_true(inherits(test_gen,"gen_tbl"))
  # we can extract the genotypes correctly
  extracted_genotypes <- test_gen %>% show_genotypes()
  expect_true(all(extracted_genotypes==test_genotypes))
  # extract them from the list directly
  expect_true(all(show_genotypes(test_gen$genotypes)==test_genotypes))
  # we can extract the loci correctly
  extracted_loci <- test_gen %>% show_loci()
  expect_identical(extracted_loci, as_tibble(test_loci))
  expect_identical( show_loci(test_gen$genotypes), as_tibble(test_loci))
  # check ploidy (it should be diploid by default)
  # not that ploidy is stored as integers
  extracted_ploidy <- test_gen %>% show_ploidy()
  expect_identical(extracted_ploidy, as.integer(rep(2,nrow(test_ind_meta))))
  expect_identical(show_ploidy(test_gen$genotypes),as.integer(rep(2,nrow(test_ind_meta))))

  # example of dropping the genotypes, leading to a change in class
  test_drop <- test_gen %>% select(-genotypes)
  expect_false(inherits(test_drop,"gen_tbl"))

})
