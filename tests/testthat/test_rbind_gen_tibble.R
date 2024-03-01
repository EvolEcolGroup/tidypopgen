#reference file
raw_path_pop_b <- system.file("extdata/pop_b.raw", package = "tidypopgen")
map_path_pop_b <- system.file("extdata/pop_b.map", package = "tidypopgen")
pop_b_gen <- read_plink_raw(file = raw_path_pop_b, map_file = map_path_pop_b, quiet = TRUE)

#target file
raw_path_pop_a <- system.file("extdata/pop_a.raw", package = "tidypopgen")
map_path_pop_a <- system.file("extdata/pop_a.map", package = "tidypopgen")
pop_a_gen <- read_plink_raw(file = raw_path_pop_a, map_file = map_path_pop_a, quiet = TRUE)

#create merge report
merged_gen <- rbind_gen_tbl(pop_a_gen,pop_b_gen, flip_strand = TRUE,
                            remove_ambiguous = TRUE, quiet = TRUE)

testthat::test_that("merge combines datasets correctly",{


  genotypes <- show_genotypes(merged_gen)

  #Genotypes before merging
  san_geno <- show_genotypes(pop_b_gen)
  sumba_geno <- show_genotypes(pop_a_gen)

  #Check San
  san_merged <- genotypes[grepl("SA", rownames(genotypes)), ]

  #All San genotypes should stay the same, as they are the reference
  testthat::expect_equal(san_merged,san_geno[,c(1,2,3,4,9,10,11,12)])

  #Check Sumba
  sumba_merged <- genotypes[grepl("GRC", rownames(genotypes)), ]

  #Check genotypes of non-swapped alleles are the same before and after merge
  testthat::expect_equal(sumba_merged[,c(1,2,4,5,6)],sumba_geno[,c(2,3,5,9,10)])

  #Check swapped alleles have genotypes swapped correctly
  testthat::expect_true(all(sumba_merged[,3] == c(2,2,0,2,1))) #rs11240777
  testthat::expect_true(all(sumba_merged[,7] == c(1,2,1,1,1))) #rs10106770
  testthat::expect_true(all(sumba_merged[,8] == c(2,1,2,1,2))) #rs11942835


  #Check ambiguous SNPs are dropped (remove_ambiguous TRUE by default)
  testthat::expect_false("rs1240719" %in% merged_gen@loc.names)
  testthat::expect_false("rs307354" %in% merged_gen@loc.names)

})

testthat::test_that("merge combines datasets correctly with file names",{
  # not yet implemented with file names
  skip(TRUE)
  #create merge report
  merged_gen_char <- em_merge(target = raw_path_pop_a,
                             ref = raw_path_pop_b,
                             target_map = map_path_pop_a,
                             ref_map = map_path_pop_b,
                             quiet = TRUE)
  # test it it identical to merged_gen
  expect_identical(merged_gen, merged_gen_char)

  # test errors
  # incorrect path
  expect_error(em_merge(target = raw_path_pop_a,
                        ref = raw_path_pop_b,
                        target_map = map_path_pop_a,
                        ref_map = "blah",
                        quiet = TRUE), "ref_map can not be")
  # missing path
  expect_error(em_merge(target = raw_path_pop_a,
                        ref = raw_path_pop_b,
                        target_map = "blah",
                        quiet = TRUE), "ref_map should be a path")


})


