testthat::test_that("select_loci subsets correctly",{
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


  # select snps with an rs
  test_gen_sub <- test_gen %>% select_loci (starts_with("rs"))
  expect_true(!any(c("x1","x2") %in% show_loci_names(test_gen_sub)))
  # subsetting by id with reordering
  test_gen_sub <- test_gen %>% select_loci (c(3,1,5))
  expect_identical(c("rs3","rs1","x1"), show_loci_names(test_gen_sub))
  # get everything
  test_gen_sub <- test_gen %>% select_loci (everything())
  expect_identical(test_gen, test_gen_sub)
  # use 2:4 range expressions
  test_gen_sub <- test_gen %>% select_loci (2:4)
  expect_identical(test_loci$name[2:4], show_loci_names(test_gen_sub))

  # now let's test some swapping
  # don't swap anything (but we use a different function for repacking that can swap)
  test_gen_sub <- test_gen %>% select_loci (c(3,1,5))
  test_gen_sub1 <- test_gen %>% select_loci (c(3,1,5),
                                             .swap_if_arg = rep(FALSE,6))
  expect_identical(test_gen_sub, test_gen_sub1)
  # swap the first SNP (which ends as the second)
  test_gen_sub2 <- test_gen %>% select_loci (c(3,1,5),
                                             .swap_if_arg = c(TRUE, rep(FALSE,5)))
  expect_true(all.equal(show_genotypes(test_gen_sub2)[,2],
                        c(1,0,0)))
  all.equal(show_loci(test_gen)[1,4:5], show_loci(test_gen_sub2)[2,5:4],
            check.attributes = FALSE)
  # swap snps based on chromosome
  test_gen_sub3 <- test_gen %>% select_loci (everything(),
                                             .swap_if_arg = show_loci(test_gen)$chromosome==2)
  expect_true(all.equal(show_genotypes(test_gen_sub3)[,5],
                        c(1,2,1)))
  all.equal(show_loci(test_gen)[5:6,4:5], show_loci(test_gen_sub3)[5:6,5:4],
            check.attributes = FALSE)
  # swap based on snp names (the "x" snps are the one on chromosome 2)
  test_gen_sub4 <- test_gen %>% select_loci (everything(),
                                             .swap_arg = starts_with("x"))
  expect_identical(test_gen_sub3,test_gen_sub4)

})
