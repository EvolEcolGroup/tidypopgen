#reference file
raw_path_pop_b <- system.file("extdata/pop_b.raw", package = "tidypopgen")
map_path_pop_b <- system.file("extdata/pop_b.map", package = "tidypopgen")
pop_b_gen <- read_plink_raw(file = raw_path_pop_b, map_file = map_path_pop_b, quiet = TRUE)

#target file
raw_path_pop_a <- system.file("extdata/pop_a.raw", package = "tidypopgen")
map_path_pop_a <- system.file("extdata/pop_a.map", package = "tidypopgen")
pop_a_gen <- read_plink_raw(file = raw_path_pop_a, map_file = map_path_pop_a, quiet = TRUE)

#create merge report
merged_gen <- rbind_gen_tbl(pop_b_gen, pop_a_gen, flip_strand = TRUE,
                            remove_ambiguous = TRUE, quiet = TRUE)

testthat::test_that("merge combines datasets correctly",{

  genotypes <- show_genotypes(merged_gen)

  #Genotypes before merging
  pop_b_geno <- show_genotypes(pop_b_gen)
  pop_a_geno <- show_genotypes(pop_a_gen)

  #Check pop_b
  pop_b_merged <- merged_gen %>% filter(population == "pop_b") %>% show_genotypes()

  #All San genotypes should stay the same, as they are the reference
  testthat::expect_equal(pop_b_merged,pop_b_geno[,c(1,2,3,4,9,10,11,12)])

  #Check pop_b
  pop_a_merged <- merged_gen %>% filter(population == "pop_a") %>% show_genotypes()

  #Check genotypes of non-swapped alleles are the same before and after merge
  testthat::expect_equal(pop_a_merged[,c(1,2,4,5,6)],pop_a_geno[,c(2,3,5,9,10)])

  #Check swapped alleles have genotypes swapped correctly
  testthat::expect_true(all(pop_a_merged[,3] == c(2,2,0,2,1))) #rs11240777
  testthat::expect_true(all(pop_a_merged[,7] == c(1,2,1,1,1))) #rs10106770
  testthat::expect_true(all(pop_a_merged[,8] == c(2,1,2,1,2))) #rs11942835


  #Check ambiguous SNPs are dropped (remove_ambiguous TRUE by default)
  testthat::expect_false("rs1240719" %in% show_loci_names(merged_gen))
  testthat::expect_false("rs307354" %in% show_loci_names(merged_gen))

})
