testthat::test_that("snpbin_list_means computes correctly",{
  test_indiv_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,2),
                          c(2,1,0,NA,0,NA),
                          c(2,2,0,0,1,NA))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=c(1,1,1,1,2,2),
                          position=c(3,5,65,343,23,456),
                          allele_ref = c("a","t","c","g","c","t"),
                          allele_alt = c("t","c", NA,"c","g","a"))
  test_gen <- gen_tibble(test_indiv_meta, test_genotypes, test_loci)

  expect_true(all(snpbin_list_means(test_gen$genotypes, use_c = FALSE)==
                    colSums(test_genotypes, na.rm=TRUE)/(c(3,3,3,2,3,1)*2)))
  expect_true(all(snpbin_list_means(test_gen$genotypes, use_c = TRUE)==
                    colSums(test_genotypes, na.rm=TRUE)/(c(3,3,3,2,3,1)*2)))
})
