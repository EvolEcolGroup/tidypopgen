

test_that("loci_alt_freq and loci_maf computes correctly",{
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

  # repeat the tests for a subset of the data
  # remove the 2nd individual and the 3rd and 5th snp
  test_genotypes <- test_genotypes[-2,c(-3,-5)]
  test_gt <- test_gt %>% filter(id!="b") %>% select_loci(c(-3,-5))
  freq <- colSums(test_genotypes, na.rm=TRUE)/(c(2,2,2,1)*2)
  expect_true(all(loci_alt_freq(test_gt$genotypes)==freq))
  # convert to minor frequencies
  freq[freq>0.5] <- 1 - freq[freq>0.5]
  expect_true(all(loci_maf(test_gt$genotypes)==freq))
})

test_that("loci_alt_freq and loci_maf on grouped tibbles",{
  test_genotypes <- rbind(c(1,1,0,1,1,0),
                          c(2,1,0,NA,0,0),
                          c(2,NA,0,0,1,1),
                          c(1,0,0,1,0,0),
                          c(1,2,0,1,2,1),
                          c(0,0,0,0,NA,1),
                          c(0,1,1,0,1,NA))
  test_indiv_meta <- data.frame (id=c("a","b","c","d","e","f","g"),
                                 population = c("pop1","pop1","pop2","pop2","pop1","pop3","pop3"))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=paste0("chr",c(1,1,1,1,2,2)),
                          position=as.integer(c(3,5,65,343,23,456)),
                          genetic_dist = as.integer(rep(0,6)),
                          allele_ref = c("A","T","C","G","C","T"),
                          allele_alt = c("T","C", NA,"C","G","A"))

  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)
  test_gt <- test_gt %>% group_by(population)
  # compute by using group map
  loci_freq_map <- test_gt %>% group_map(.f=~loci_alt_freq(.x))
  # use fast cpp code (limit cores to 2)
  loci_freq_grp <- test_gt %>% loci_alt_freq(n_cores=2)
  all.equal(loci_freq_map, loci_freq_grp)
  # and now for maf
  # compute by using group map
  loci_maf_map <- test_gt %>% group_map(.f=~loci_maf(.x))
  # use fast cpp code (limit cores to 2)
  loci_maf_grp <- test_gt %>% loci_maf_freq(n_cores=2)
  all.equal(loci_maf_map, loci_maf_grp)
})
