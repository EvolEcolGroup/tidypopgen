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
test_that("save and load gt",{
  expect_true(inherits(test_gt,"gen_tbl"))
  # now save the tibble
  all_file_names <- gt_save(test_gt)
  # check that the new file exists
  expect_true(file.exists(all_file_names[1]))
  new_test_gt <- gt_load(all_file_names[1])

  # check that we preserved the genotypes
  expect_true(all(show_genotypes(new_test_gt$genotypes)==test_genotypes))
  # check that we preserved the loci
  expect_identical( show_loci(new_test_gt$genotypes) %>% select(-big_index), as_tibble(test_loci))

})


