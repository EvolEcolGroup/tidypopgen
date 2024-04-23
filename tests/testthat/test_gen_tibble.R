# create file
test_indiv_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
test_genotypes <- rbind(c(1,1,0,1,1,0),
                        c(2,1,0,0,0,0),
                        c(2,2,0,0,1,1))
test_loci <- data.frame(name=paste0("rs",1:6),
                        chromosome=paste0("chr",c(1,1,1,1,2,2)),
                        position=as.integer(c(3,5,65,343,23,456)),
                        genetic_dist = as.integer(rep(0,6)),
                        allele_ref = c("A","T","C","G","C","T"),
                        allele_alt = c("T","C", NA,"C","G","A"))

test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)

# this also tests show_genotypes and show_loci
test_that("create gen_tibble from bed",{
  expect_true(inherits(test_gt,"gen_tbl"))
  # we can extract the genotypes correctly
  extracted_genotypes <- test_gt %>% show_genotypes()
  expect_true(all(extracted_genotypes==test_genotypes))
  # extract them from the list directly
  expect_true(all(show_genotypes(test_gt$genotypes)==test_genotypes))
  # we can extract the loci correctly
  extracted_loci <- test_gt %>% show_loci()
  # remove the index in the big file
  expect_identical( show_loci(test_gt$genotypes) %>% select(-big_index), as_tibble(test_loci))

  # example of dropping the genotypes, leading to a change in class
  test_drop <- test_gt %>% select(-genotypes)
  expect_false(inherits(test_drop,"gen_tbl"))

})

# now create it directly from the dfs
test_that("create gen_tibble from dfs",{
  test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = test_indiv_meta,
             loci = test_loci, quiet = TRUE)
  # because of the different backing file info, we cannot use identical on the whole object
  expect_true(identical(show_genotypes(test_gt), show_genotypes(test_dfs_gt)))
  expect_true(identical(show_loci(test_gt), show_loci(test_dfs_gt)))
  expect_true(identical(test_gt %>% select(-genotypes),
                        test_dfs_gt %>% select(-genotypes)))
})

test_genotypes_c <- rbind(c("1","1","0","1","1","0"),
                          c("2","1","0","0","0","0"),
                          c("2","2","0","0","1","1"))


test_that("gen_tibble does not accept character matrix",{
  expect_error(test_dfs_gt <- gen_tibble(test_genotypes_c, indiv_meta = test_indiv_meta,
                                         loci = test_loci, quiet = TRUE),"'x' is not a matrix of integers")
})

test_that("gen_tibble catches invalid alleles",{
  test_loci_wrong <- test_loci
  test_loci_wrong$allele_alt[1] <- "N"
  expect_error(test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = test_indiv_meta,
                                           loci = test_loci_wrong, quiet = TRUE),"valid alleles are")
  # now add N to the valid alleles
  test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = test_indiv_meta,
                                         loci = test_loci_wrong,
                                         valid_alleles = c("A","C","T","G","N"),
                            quiet = TRUE)
  expect_true("N" %in% show_loci(test_dfs_gt)$allele_alt)
  # but if we add to missing values it shoudl be turned into a zero
  test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = test_indiv_meta,
                            loci = test_loci_wrong,
                            missing_alleles = c("0",".","N"),
                            quiet = TRUE)
  expect_false("N" %in% show_loci(test_dfs_gt)$allele_alt)
  expect_true(is.na(show_loci(test_dfs_gt)$allele_alt[1]))
  # and finally throw an error if we try to use 0 as a missing value
  expect_error(test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = test_indiv_meta,
                            loci = test_loci_wrong,
                            valid_alleles = c("A","C","T","G","0"),
                            quiet = TRUE), "can not be a valid allele")



})


