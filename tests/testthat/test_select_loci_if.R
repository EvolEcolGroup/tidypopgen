testthat::test_that(".genotypes_means computes correctly",{
  test_ind_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,2),
                          c(2,1,0,NA,0,NA),
                          c(2,2,0,0,1,NA))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=c(1,1,1,1,2,2),
                          position=c(3,5,65,343,23,456),
                          allele_ref = c("a","t","c","g","c","t"),
                          allele_alt = c("t","c", NA,"c","g","a"))
  test_gen <- gen_tibble(test_ind_meta, test_genotypes, test_loci)

  # raw frequencies
  freq <- test_gen %>% loci_freq()
  # remove monomorphic
  test_gen_sub <- test_gen %>% select_loci_if(loci_freq(genotypes)!=0)
  expect_true(!any(c("rs3","rs6") %in% show_loci(test_gen_sub)$name))
  # same subsetting with .data
  expect_identical(test_gen_sub,
                   test_gen %>% select_loci_if(loci_freq(.data$genotypes)!=0))
  # same manually
  criterion <- c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)
  expect_identical(test_gen_sub,
                   test_gen %>% select_loci_if(criterion))
  #TODO test errors
})
