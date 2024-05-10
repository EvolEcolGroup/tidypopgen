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
test_that("create gen_tibble from dfs",{
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

test_that("gen_tibble from a bed file",{
  bed_path <- system.file("extdata/pop_a.bed", package = "tidypopgen")
  pop_a_gt <- gen_tibble(bed_path, quiet=TRUE, backingfile = tempfile())
  # now read the dosages created by plink when saving in raw format
  raw_file_pop_a <- read.table(system.file("extdata/pop_a.raw", package = "tidypopgen"), header= TRUE)
  mat <- as.matrix(raw_file_pop_a[,7:ncol(raw_file_pop_a)])
  mat <- unname(mat)
  expect_true(all.equal(mat,show_genotypes(pop_a_gt)))
  # now read in the ped file
  ped_path <- system.file("extdata/pop_a.ped", package = "tidypopgen")
  pop_a_ped_gt <- gen_tibble(ped_path, quiet=TRUE,backingfile = tempfile())
  all.equal(show_genotypes(pop_a_gt),show_genotypes(pop_a_ped_gt))

})
