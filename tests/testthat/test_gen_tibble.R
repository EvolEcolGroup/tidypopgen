source("make_test_bed")
# this also tests show_genotypes and show_loci
testthat::test_that("gen_tibble stores data correctly",{
  # create file
  test_ind_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,0),
                          c(2,1,0,0,0,0),
                          c(2,2,0,0,1,1))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=paste0("chr",c(1,1,1,1,2,2)),
                          position=as.integer(c(3,5,65,343,23,456)),
                          genetic_dist = as.integer(rep(0,6)),
                          allele_ref = c("a","t","c","g","c","t"),
                          allele_alt = c("t","c", NA,"c","g","a"))

  file_plink <- tempfile('test_data_')
  make_test_bed(test_genotypes, test_loci, test_ind_meta, file_plink)
  file_plink<-paste0(file_plink,".bed")
  # convert bed to bigsnp
  path_rds <- bigsnpr::snp_readBed(file_plink, backingfile = tempfile("test_bigfile_"))
  # convert to gen_tibble
  test_gt <- gen_tibble(path_rds)
  expect_true(inherits(test_gt,"gen_tbl"))
  # we can extract the genotypes correctly
  extracted_genotypes <- test_gt %>% show_genotypes()
  expect_true(all(extracted_genotypes==test_genotypes))
  # extract them from the list directly
  expect_true(all(show_genotypes(test_gt$genotypes)==test_genotypes))
  # we can extract the loci correctly
  extracted_loci <- test_gt %>% show_loci()
  expect_identical(extracted_loci, as_tibble(test_loci))
  # remove the index in the big file
  expect_identical( show_loci(test_gt$genotypes) %>% select(-big_index), as_tibble(test_loci))

  # example of dropping the genotypes, leading to a change in class
  test_drop <- test_gt %>% select(-genotypes)
  expect_false(inherits(test_drop,"gen_tbl"))

})
