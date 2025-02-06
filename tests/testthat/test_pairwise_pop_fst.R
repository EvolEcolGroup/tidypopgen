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
                        genetic_dist = as.double(rep(0,6)),
                        allele_ref = c("A","T","C","G","C","T"),
                        allele_alt = c("T","C", NA,"C","G","A"))

test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)


test_that("pairwise_pop_fst compute correctly",{
  test_gt <- test_gt %>% dplyr::group_by(population)
  test_hier <- gt_as_hierfstat(test_gt)
  # compare results against hierfstat for Nei87 (Nei86 does not correct for Ho
  # when computing Ht, so it gives a different result)
  nei_gt <- test_gt %>% pairwise_pop_fst(method="Nei87")

  nei_hier <- hierfstat::pairwise.neifst(test_hier)
  # hiefstat values are rounded to 4 dp
  expect_true(all.equal(tidy_dist_matrix(nei_hier)$value, round(nei_gt$value,4)))

  #pair_fst_locus <- test_gt %>% pairwise_pop_fst(by_locus = TRUE)

})


test_that("pairwise_pop_fst weir_cockerham",{
  test_genotypes <- rbind(c(1,1,0,1,1,0),
                          c(2,1,0,NA,0,0),
                          c(2,NA,0,0,1,1),
                          c(1,0,0,1,0,0),
                          c(1,2,0,1,2,1),
                          c(0,0,0,0,NA,1),
                          c(0,1,1,0,1,NA))
  test_indiv_meta <- data.frame (id=c("a","b","c","d","e","f","g"),
                                 population = c("pop1","pop1","pop1","pop2","pop2","pop2","pop2"))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=paste0("chr",c(1,1,1,1,2,2)),
                          position=as.integer(c(3,5,65,343,23,456)),
                          genetic_dist = as.double(rep(0,6)),
                          allele_ref = c("A","T","C","G","C","T"),
                          allele_alt = c("T","C", NA,"C","G","A"))
  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)

  test_gt <- test_gt %>% dplyr::group_by(population)
  wc_tidypopgen <- test_gt %>% pairwise_pop_fst(method="WC84")

  # compared to hierfstat
  test_heir <- gt_as_hierfstat(test_gt)
  wc_hierfstat <- hierfstat::pairwise.WCfst(test_heir)
  expect_equal(wc_tidypopgen$value, wc_hierfstat[1,2])

  # compared to scikit-allel version 1.3.13
  # See create_scikit-allel_test_data for script

  wc_scikit <- as.numeric(readLines(system.file("/extdata/fst_scikit-allel/fst_wc.txt", package = "tidypopgen")))
  expect_equal(wc_tidypopgen$value, wc_scikit)

  ########### test with a monomorphic loci
  test_genotypes <- rbind(c(2,1,0,1,1,0),
                          c(2,1,0,NA,0,0),
                          c(2,NA,0,0,1,1),
                          c(2,0,0,1,0,0),
                          c(2,2,0,1,2,1),
                          c(2,0,0,0,NA,1),
                          c(2,1,1,0,1,NA))
  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)

  test_gt <- test_gt %>% dplyr::group_by(population)
  wc_tidypopgen_mono <- test_gt %>% pairwise_pop_fst(method="WC84")

  # read in output
  wc_scikit_mono <- as.numeric(readLines(system.file("/extdata/fst_scikit-allel/fst_wc_monomorphic.txt", package = "tidypopgen")))
  expect_equal(wc_tidypopgen_mono$value, wc_scikit_mono)
})


test_that("pairwise_pop_fst hudson",{
  test_genotypes <- rbind(c(1,1,0,1,1,0),
                          c(2,1,0,NA,0,0),
                          c(2,NA,0,0,1,1),
                          c(1,0,0,1,0,0),
                          c(1,2,0,1,2,1),
                          c(0,0,0,0,NA,1),
                          c(0,1,1,0,1,NA))
  test_indiv_meta <- data.frame (id=c("a","b","c","d","e","f","g"),
                                 population = c("pop1","pop1","pop1","pop2","pop2","pop2","pop2"))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=paste0("chr",c(1,1,1,1,2,2)),
                          position=as.integer(c(3,5,65,343,23,456)),
                          genetic_dist = as.double(rep(0,6)),
                          allele_ref = c("A","T","C","G","C","T"),
                          allele_alt = c("T","C", NA,"C","G","A"))
  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)

  test_gt <- test_gt %>% dplyr::group_by(population)
  hudson_tidypopgen <- test_gt %>% pairwise_pop_fst(method="Hudson")

  # compared to scikit-allel version 1.3.13
  # See create_scikit-allel_test_data for script

  # read in output
  hudson_scikit <- as.numeric(readLines(system.file("/extdata/fst_scikit-allel/fst_hudson.txt", package = "tidypopgen")))
  expect_equal(hudson_tidypopgen$value, hudson_scikit)

  ########### test with a monomorphic loci
  test_genotypes <- rbind(c(2,1,0,1,1,0),
                          c(2,1,0,NA,0,0),
                          c(2,NA,0,0,1,1),
                          c(2,0,0,1,0,0),
                          c(2,2,0,1,2,1),
                          c(2,0,0,0,NA,1),
                          c(2,1,1,0,1,NA))
  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)

  test_gt <- test_gt %>% dplyr::group_by(population)
  hudson_tidypopgen_mono <- test_gt %>% pairwise_pop_fst(method="Hudson")

  # read in output
  hudson_scikit_mono <- as.numeric(readLines(system.file("/extdata/fst_scikit-allel/fst_hudson_monomorphic.txt", package = "tidypopgen")))
  expect_equal(hudson_tidypopgen_mono$value, hudson_scikit_mono)

})


