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

test_that("warning when no snps overlap",{

  #Create two datasets

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


  test_indiv_meta2 <- data.frame (id=c("a","b","c"),
                                  population = c("pop1","pop1","pop2"))
  test_genotypes2 <- rbind(c(1,1,0,1,1,2),
                           c(2,1,0,NA,0,NA),
                           c(2,2,0,0,1,NA))
  test_loci2 <- data.frame(name=paste0("rs",7:12),
                           chromosome=c(1,1,1,1,2,2),
                           position=c(3,5,65,343,23,456),
                           genetic_dist = as.integer(rep(0,6)),
                           allele_ref = c("A","T","C","G","C","T"),
                           allele_alt = c("T","C", NA,"C","G","A"))
  test_gt2 <- gen_tibble(x = test_genotypes2, loci = test_loci2, indiv_meta = test_indiv_meta2, quiet = TRUE)

  #merge
  expect_error(rbind.gen_tbl(test_gt, test_gt2, flip_strand = TRUE,
                             quiet = TRUE,
                             backingfile = tempfile()), "There are no overlapping loci")

  #Now try with the same snp ID but different CHR

  #These should merge, with CHR and POS in merged GT being equal to test_gt2

  test_indiv_meta2 <- data.frame (id=c("a","b","c"),
                                  population = c("pop1","pop1","pop2"))
  test_genotypes2 <- rbind(c(1,1,0,1,1,2),
                           c(2,1,0,NA,0,NA),
                           c(2,2,0,0,1,NA))
  test_loci2 <- data.frame(name=paste0("rs",1:6),
                           chromosome=c(3,3,3,3,4,4),
                           position=c(3,5,65,343,23,456),
                           genetic_dist = as.integer(rep(0,6)),
                           allele_ref = c("A","T","C","G","C","T"),
                           allele_alt = c("T","C", NA,"C","G","A"))
  test_gt2 <- gen_tibble(x = test_genotypes2, loci = test_loci2, indiv_meta = test_indiv_meta2, quiet = TRUE)

  report <- rbind_dry_run(test_gt, test_gt2, flip_strand = TRUE)

  #merge
  merged_gt <- rbind.gen_tbl(test_gt, test_gt2, flip_strand = TRUE,
                             quiet = TRUE,
                             backingfile = tempfile())

  #only two remain - ambiguous are removed but message suggests they aren't?
  #"( 4 are ambiguous, of which 0 were removed)"
  expect_true(all(show_loci(merged_gt)$"chromosome" == c(1,1)))
  expect_true(all(show_loci(merged_gt)$"position" == c(5,65)))

})



