raw_path_pop_b <- system.file("extdata/pop_b.bed", package = "tidypopgen")
bigsnp_path_b <- bigsnpr::snp_readBed(raw_path_pop_b, backingfile = tempfile("test_b_"))
pop_b_gt <- gen_tibble(bigsnp_path_b, quiet=TRUE)
#target file
raw_path_pop_a <- system.file("extdata/pop_a.bed", package = "tidypopgen")
bigsnp_path_a <- bigsnpr::snp_readBed(raw_path_pop_a, backingfile = tempfile("test_a_"))
pop_a_gt <- gen_tibble(bigsnp_path_a, quiet=TRUE)
# #create merge
 merged_gen <- rbind.gen_tbl(pop_b_gt, pop_a_gt, flip_strand = TRUE,
                             quiet = TRUE,
                             backingfile = tempfile())

test_that("merge combines datasets correctly",{

  genotypes <- show_genotypes(merged_gen)

  #Genotypes before merging
  pop_b_geno <- show_genotypes(pop_b_gt)
  pop_a_geno <- show_genotypes(pop_a_gt)

  #Check pop_b
  pop_b_merged <- merged_gen %>% filter(population == "pop_b") %>% show_genotypes()

  #All pop_b genotypes should stay the same, as they are the reference
  testthat::expect_equal(pop_b_merged,pop_b_geno[,c(1,2,3,4,5,6,7,13,14,15,16,17)])

  #Check pop_b
  pop_a_merged <- merged_gen %>% filter(population == "pop_a") %>% show_genotypes()

  #Check genotypes of non-swapped alleles are the same before and after merge
  testthat::expect_equal(pop_a_merged[,c(1,2,3)],pop_a_geno[,c(2,3,4)])

  #Check swapped alleles have genotypes swapped correctly
  testthat::expect_true(all(pop_a_merged[,4] == c(2,2,0,2,1))) #rs11240777
  testthat::expect_true(all(pop_a_merged[,10] == c(1,2,1,1,1))) #rs10106770
  testthat::expect_true(all(pop_a_merged[,11] == c(2,1,2,1,2))) #rs11942835


  #Check ambiguous SNPs are dropped (remove_ambiguous TRUE by default)
  testthat::expect_false("rs1240719" %in% loci_names(merged_gen))
  testthat::expect_false("rs307354" %in% loci_names(merged_gen))

})
