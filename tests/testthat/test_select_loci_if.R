testthat::test_that("select_loci_if subsets correctly",{
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


  # now let's test some swapping
  # don't swap anything (but we use a different function for repacking that can swap)
  test_gen_sub <- test_gen %>% select_loci_if (criterion)
  test_gen_sub1 <- test_gen %>% select_loci_if (criterion,
                                             .swap_if_arg = rep(FALSE,6))
  expect_identical(test_gen_sub, test_gen_sub1)
  # swap the first SNP
  test_gen_sub2 <- test_gen %>% select_loci_if (criterion,
                                             .swap_if_arg = c(TRUE, rep(FALSE,5)))
  expect_true(all.equal(show_genotypes(test_gen_sub2)[,1],
                        c(1,0,0)))
  # check the alleles
  all.equal(show_loci(test_gen)[1,4:5], show_loci(test_gen_sub2)[1,5:4],
            check.attributes = FALSE)
  # swap snps based on chromosome
  test_gen_sub3 <- test_gen %>% select_loci_if (rep(TRUE, 6),
                                             .swap_if_arg = show_loci(test_gen)$chromosome==2)
  expect_true(all.equal(show_genotypes(test_gen_sub3)[,5],
                        c(1,2,1)))
  all.equal(show_loci(test_gen)[5:6,4:5], show_loci(test_gen_sub3)[5:6,5:4],
            check.attributes = FALSE)

})
