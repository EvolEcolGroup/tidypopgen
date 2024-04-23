test_that("snpbin_list_means computes correctly",{
  test_indiv_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,2),
                          c(2,1,0,NA,0,NA),
                          c(2,2,0,0,1,NA))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=c(1,1,1,1,2,2),
                          position=c(3,5,65,343,23,456),
                          genetic_dist = as.integer(rep(0,6)),
                          allele_ref = c("A","T","C","G","C","T"),
                          allele_alt = c("T","C", NA,"C","G","A"))

  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)


  # raw frequencies
  freq <- colSums(test_genotypes, na.rm=TRUE)/(c(3,3,3,2,3,1)*2)
  expect_true(all(loci_alt_freq(test_gt$genotypes)==freq))
  # convert to minor frequencies
  freq[freq>0.5] <- 1 - freq[freq>0.5]
  expect_true(all(loci_maf(test_gt$genotypes)==freq))

})
