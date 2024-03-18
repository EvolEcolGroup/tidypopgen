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
bed_path <- gt_write_bed_from_dfs(genotypes = test_genotypes,
                                  loci = test_loci,
                                  indiv_meta = test_indiv_meta,
                                  path_out = tempfile('test_data_'))
test_gt <- gen_tibble(bed_path, quiet = TRUE)

# we now replace NA with 0 for the test_loci
test_loci[is.na(test_loci)]<-"0"

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

