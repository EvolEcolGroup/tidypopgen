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


test_that("if order of loci is changed, order of genotypes also changes",{

  pop_b <- gen_tibble(system.file("extdata/pop_b.bed", package="tidypopgen"),backingfile = tempfile(), quiet = TRUE)
  #original genotypes
  pop_b_gen <- show_genotypes(pop_b)

  #now scramble the loci
  set.seed(123)
  random_order <- sample(1:17)
  show_loci(pop_b) <- pop_b %>% select_loci(all_of(random_order)) %>% show_loci()

  #reorder the original genotypes according to 'random_order'
  pop_b_gen_reordered <- pop_b_gen[,random_order]

  #check that genotypes are now reordered according to random order
  expect_equal(pop_b_gen_reordered, show_genotypes(pop_b))


})

test_that("gen_tibble does not accept character matrix",{
  test_genotypes_c <- rbind(c("1","1","0","1","1","0"),
                            c("2","1","0","0","0","0"),
                            c("2","2","0","0","1","1"))
  expect_error(test_dfs_gt <- gen_tibble(test_genotypes_c, indiv_meta = test_indiv_meta,
                                         loci = test_loci, quiet = TRUE),"'x' is not a matrix of integers")
})

test_that("gen_tibble wrong filetype error",{
  expect_error(test_dfs_gt <- gen_tibble(system.file("extdata/related/test_king.kin0", package = "tidypopgen")),
               "file_path should be pointing ")
})

test_that("gen_tibble loci is dataframe or tbl",{

  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=paste0("chr",c(1,1,1,1,2,2)),
                          position=as.integer(c(3,5,65,343,23,456)),
                          genetic_dist = as.integer(rep(0,6)),
                          allele_ref = c("A","T","C","G","C","T"),
                          allele_alt = c("T","C", NA,"C","G","A"))
  wrong_loci_matrix <- as.matrix(test_loci)

  expect_error(test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = test_indiv_meta,
                                         loci = wrong_loci_matrix, quiet = TRUE),"loci must be one of data.frame or tbl")
})

test_that("gen_tibble required id and population",{
  wrong_indiv_meta <- data.frame (x =c("a","b","c"),
                                  y = c("pop1","pop1","pop2"))
  expect_error(test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = wrong_indiv_meta,
                                         loci = test_loci, quiet = TRUE),"ind_meta does not include the compulsory columns")
})

test_that("gen_tibble indiv_meta is list, dataframe, or tbl",{
  wrong_indiv_meta <- data.frame (id=c("a","b","c"),
                                  population = c("pop1","pop1","pop2"))
  wrong_indiv_meta_matrix <- as.matrix(wrong_indiv_meta)

  expect_error(test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = wrong_indiv_meta_matrix,
                                         loci = test_loci, quiet = TRUE),"indiv_meta must be one of data.frame, tbl, or list")
})

test_that("gen_tibble identifies wrong dimensions in genotypes",{
  wrong_genotypes <- rbind(c(1,1,0,1,1,0),
                          c(2,1,0,0,0,0))
  expect_error(test_dfs_gt <- gen_tibble(wrong_genotypes, indiv_meta = test_indiv_meta,
                                         loci = test_loci, quiet = TRUE),
               "there is a mismatch between the number of loci in the genotype table x and in the loci table")

})

test_that("gen_tibble identifies wrong loci table columns",{
  wrong_loci <- data.frame(a=paste0("rs",1:6),
                          b=paste0("chr",c(1,1,1,1,2,2)),
                          c=as.integer(c(3,5,65,343,23,456)),
                          d = as.integer(rep(0,6)),
                          e = c("A","T","C","G","C","T"),
                          f = c("T","C", NA,"C","G","A"))
  expect_error(test_dfs_gt <- gen_tibble(test_genotypes, indiv_meta = test_indiv_meta,
                                         loci = wrong_loci, quiet = TRUE),
               "loci does not include the compulsory columns")
})


test_that("gen_tibble from files",{
  ########################
  # PLINK BED files
  ########################
  bed_path <- system.file("extdata/pop_a.bed", package = "tidypopgen")
  pop_a_gt <- gen_tibble(bed_path, quiet=TRUE, backingfile = tempfile())
  # now read the dosages created by plink when saving in raw format
  raw_file_pop_a <- read.table(system.file("extdata/pop_a.raw", package = "tidypopgen"), header= TRUE)
  mat <- as.matrix(raw_file_pop_a[,7:ncol(raw_file_pop_a)])
  mat <- unname(mat)
  expect_true(all.equal(mat,show_genotypes(pop_a_gt)))
  ########################
  # PLINK PED files
  ########################
  ped_path <- system.file("extdata/pop_a.ped", package = "tidypopgen")
  pop_a_ped_gt <- gen_tibble(ped_path, quiet=TRUE,backingfile = tempfile())
  # because ref and alt are defined based on which occurs first in a ped, some alleles will be swapped
  equal_geno <- show_genotypes(pop_a_gt)==show_genotypes(pop_a_ped_gt)
  not_equal <- which(!apply(equal_geno,2,all))
  # check that the alleles for loci that are mismatched are indeed swapped
  expect_true(all(show_loci(pop_a_gt)$allele_alt[not_equal] == show_loci(pop_a_ped_gt)$allele_ref[not_equal]))
  # check that the mismatches are all in the homozygotes
  expect_true(all(abs(show_genotypes(pop_a_gt)[, not_equal]-show_genotypes(pop_a_ped_gt)[, not_equal]) %in% c(0,2)))
  ########################
  # PLINK VCF files
  ########################
  vcf_path <- system.file("extdata/pop_a.vcf", package = "tidypopgen")
  pop_a_vcf_gt <- gen_tibble(vcf_path, quiet=TRUE,backingfile = tempfile())
  expect_true(all.equal(show_genotypes(pop_a_gt),show_genotypes(pop_a_vcf_gt)))
  # reload it in chunks
  pop_a_vcf_gt2 <- gen_tibble(vcf_path, quiet=TRUE,backingfile = tempfile(),
                              chunk_size = 2)
  expect_true(all.equal(show_genotypes(pop_a_vcf_gt2),show_genotypes(pop_a_vcf_gt)))
  expect_true(all.equal(show_loci(pop_a_vcf_gt2),show_loci(pop_a_vcf_gt)))

  # @TODO we should add similar tests for pop b, which has missing data

  # check our cpp parser
  pop_a_vcf_fast_gt <- gen_tibble(vcf_path, quiet=TRUE,backingfile = tempfile(), parser="cpp")
  expect_true(all.equal(show_genotypes(pop_a_gt),show_genotypes(pop_a_vcf_fast_gt)))
  # check loci table against the vcfR parser
  expect_true(all.equal(show_loci(pop_a_vcf_gt), show_loci(pop_a_vcf_fast_gt)))
  # reload it in chunks
  pop_a_vcf_fast_gt2 <- gen_tibble(vcf_path, quiet=TRUE, backingfile = tempfile(),
                              chunk_size = 2, parser="cpp")
  expect_true(all.equal(show_genotypes(pop_a_vcf_fast_gt2),show_genotypes(pop_a_vcf_fast_gt)))
  expect_true(all.equal(show_loci(pop_a_vcf_gt), show_loci(pop_a_vcf_fast_gt)))

})

test_that("gen_tibble from files with missingness",{
  ########################
  # PLINK BED files
  ########################
  bed_path <- system.file("extdata/pop_b.bed", package = "tidypopgen")
  pop_b_gt <- gen_tibble(bed_path, quiet=TRUE, backingfile = tempfile())
  # now read the dosages created by plink when saving in raw format
  raw_file_pop_b <- read.table(system.file("extdata/pop_b.raw", package = "tidypopgen"), header= TRUE)
  mat <- as.matrix(raw_file_pop_b[,7:ncol(raw_file_pop_b)])
  mat <- unname(mat)
  expect_true(all.equal(mat,show_genotypes(pop_b_gt)))
  ########################
  # PLINK PED files
  ########################
  ped_path <- system.file("extdata/pop_b.ped", package = "tidypopgen")
  pop_b_ped_gt <- gen_tibble(ped_path, quiet=TRUE,backingfile = tempfile())
  # because ref and alt are defined based on which occurs first in a ped, some alleles will be swapped
  equal_geno <- show_genotypes(pop_b_gt)==show_genotypes(pop_b_ped_gt)
  not_equal <- which(!apply(equal_geno,2,all))
  # check that the alleles for loci that are mismatched are indeed swapped
  expect_true(all(show_loci(pop_b_gt)$allele_alt[not_equal] == show_loci(pop_b_ped_gt)$allele_ref[not_equal]))
  # check that the mismatches are all in the homozygotes
  expect_true(all(abs(show_genotypes(pop_b_gt)[, not_equal]-show_genotypes(pop_b_ped_gt)[, not_equal]) %in% c(0,2)))
  ########################
  # PLINK VCF files
  ########################
  vcf_path <- system.file("extdata/pop_b.vcf", package = "tidypopgen")
  pop_b_vcf_gt <- gen_tibble(vcf_path, quiet=TRUE,backingfile = tempfile())
  expect_true(all.equal(show_genotypes(pop_b_gt),show_genotypes(pop_b_vcf_gt)))
  # reload it in chunks
  pop_b_vcf_gt2 <- gen_tibble(vcf_path, quiet=TRUE,backingfile = tempfile(),
                              chunk_size = 2)
  expect_true(all.equal(show_genotypes(pop_b_vcf_gt2),show_genotypes(pop_b_vcf_gt)))
  expect_true(all.equal(show_loci(pop_b_vcf_gt2),show_loci(pop_b_vcf_gt)))

})

test_that("gentibble with packedancestry",{
  geno_path <- system.file("extdata/pop_a.geno", package = "tidypopgen")
  pop_a_gt <- gen_tibble(geno_path, quiet=TRUE, backingfile = tempfile(), valid_alleles = c("A","G","C","T"))

  #dosages in packedancestry are the opposite to .raw
  raw_file_pop_a <- read.table(system.file("extdata/pop_a.raw", package = "tidypopgen"), header= TRUE)
  mat <- as.matrix(raw_file_pop_a[,7:ncol(raw_file_pop_a)])
  mat <- unname(mat)

  zero_positions <- mat == 0
  two_positions <- mat == 2

  #swap around the dosages in .raw matrix to check correspondence with packedancestry genotypes
  mat[zero_positions] <- -1
  mat[two_positions] <- 0
  mat[mat == -1] <- 2
  expect_true(all.equal(mat,show_genotypes(pop_a_gt)))

  #read in ped to check loci
  ped_path <- system.file("extdata/pop_a.ped", package = "tidypopgen")
  pop_a_gt_ped <- gen_tibble(ped_path, quiet=TRUE, backingfile = tempfile(), valid_alleles = c("A","G","C","T"))

  #allele_alt and allele_ref are also swapped in packedancestry

  #snps 6 and 14 both have their genotypes and loci switched
  #rs1110052 and rs10106770
  #for these two, the gen_tibble from ped and gen_tibble from geno should therefore match

  keep <- which(show_loci(pop_a_gt_ped)$name %in% c("rs1110052","rs10106770"))
  pop_a_gt_ped_swap <- pop_a_gt_ped %>% select_loci(all_of(keep))
  keep <- which(show_loci(pop_a_gt)$name %in% c("rs1110052","rs10106770"))
  pop_a_gt_swap <- pop_a_gt %>% select_loci(all_of(keep))

  #check loci and genotypes
  expect_equal(show_loci(pop_a_gt_ped_swap)$allele_ref, show_loci(pop_a_gt_swap)$allele_ref)
  expect_equal(show_loci(pop_a_gt_ped_swap)$allele_alt, show_loci(pop_a_gt_swap)$allele_alt)
  expect_equal(show_genotypes(pop_a_gt_ped_swap), show_genotypes(pop_a_gt_swap))


  #for the rest, allele ref and allele alt should match
  keep <- which(!show_loci(pop_a_gt_ped)$name %in% c("rs1110052","rs10106770"))
  pop_a_gt_ped_snp <- pop_a_gt_ped %>% select_loci(all_of(keep))
  keep <- which(!show_loci(pop_a_gt)$name %in% c("rs1110052","rs10106770"))
  pop_a_gt_snp <- pop_a_gt %>% select_loci(all_of(keep))

  expect_equal(show_loci(pop_a_gt_ped_snp)$allele_ref, show_loci(pop_a_gt_snp)$allele_alt)
  expect_equal(show_loci(pop_a_gt_ped_snp)$allele_alt, show_loci(pop_a_gt_snp)$allele_ref)

  #0 and 2 dosages will be switched
  ped_gtypes <- show_genotypes(pop_a_gt_ped)
  ped_gtypes <- unname(ped_gtypes)
  zero_positions <- ped_gtypes == 0
  two_positions <- ped_gtypes == 2

  #swap around the dosages in .ped_gtypes matrix to check correspondence with packedancestry genotypes
  ped_gtypes[zero_positions] <- -1
  ped_gtypes[two_positions] <- 0
  ped_gtypes[ped_gtypes == -1] <- 2

  #remove the snps above
  ped_gtypes <- ped_gtypes[,-c(6,14)]
  expect_equal(ped_gtypes,show_genotypes(pop_a_gt_snp))

})

test_that("gentibble with packedancestry and missingness",{
  geno_path <- system.file("extdata/pop_b.geno", package = "tidypopgen")
  pop_b_gt <- gen_tibble(geno_path, quiet=TRUE, backingfile = tempfile(),
                         valid_alleles = c("A","G","C","T"))

  #dosages in packedancestry are the opposite to .raw
  raw_file_pop_b <- read.table(system.file("extdata/pop_b.raw", package = "tidypopgen"), header= TRUE)
  mat <- as.matrix(raw_file_pop_b[,7:ncol(raw_file_pop_b)])
  mat <- unname(mat)

  zero_positions <- mat == 0
  two_positions <- mat == 2

  #swap around the dosages in .raw matrix to check correspondence with packedancestry genotypes
  mat[zero_positions] <- -1
  mat[two_positions] <- 0
  mat[mat == -1] <- 2
  expect_true(all.equal(mat,show_genotypes(pop_b_gt)))

  #read in ped to check loci
  ped_path <- system.file("extdata/pop_b.ped", package = "tidypopgen")
  pop_b_gt_ped <- gen_tibble(ped_path, quiet=TRUE, backingfile = tempfile(),
                             valid_alleles = c("A","G","C","T"))

  #allele_alt and allele_ref are also swapped in packedancestry
  #snps 15 and 16 both have their genotypes and loci switched
  #rs10106770 and rs11942835
  #for these two, the gen_tibble from ped and gen_tibble from geno should therefore match

  keep <- which(show_loci(pop_b_gt_ped)$name %in% c("rs11942835","rs10106770"))
  pop_b_gt_ped_swap <- pop_b_gt_ped %>% select_loci(all_of(keep))
  keep <- which(show_loci(pop_b_gt)$name %in% c("rs11942835","rs10106770"))
  pop_b_gt_swap <- pop_b_gt %>% select_loci(all_of(keep))

  #check loci and genotypes
  expect_equal(show_loci(pop_b_gt_ped_swap)$allele_ref, show_loci(pop_b_gt_swap)$allele_ref)
  expect_equal(show_loci(pop_b_gt_ped_swap)$allele_alt, show_loci(pop_b_gt_swap)$allele_alt)
  expect_equal(show_genotypes(pop_b_gt_ped_swap), show_genotypes(pop_b_gt_swap))

  #for the rest, allele ref and allele alt should match
  keep <- which(!show_loci(pop_b_gt_ped)$name %in% c("rs11942835","rs10106770"))
  pop_b_gt_ped_snp <- pop_b_gt_ped %>% select_loci(all_of(keep))
  keep <- which(!show_loci(pop_b_gt)$name %in% c("rs11942835","rs10106770"))
  pop_b_gt_snp <- pop_b_gt %>% select_loci(all_of(keep))

  #check loci and genotypes
  expect_equal(show_loci(pop_b_gt_ped_snp)$allele_ref, show_loci(pop_b_gt_snp)$allele_alt)
  expect_equal(show_loci(pop_b_gt_ped_snp)$allele_alt, show_loci(pop_b_gt_snp)$allele_ref)

  #and genotypes should be the opposite

  #0 and 2 dosages will be switched
  ped_gtypes <- show_genotypes(pop_b_gt_ped)
  ped_gtypes <- unname(ped_gtypes)
  zero_positions <- ped_gtypes == 0
  two_positions <- ped_gtypes == 2

  #swap around the dosages in .ped_gtypes matrix to check correspondence with packedancestry genotypes
  ped_gtypes[zero_positions] <- -1
  ped_gtypes[two_positions] <- 0
  ped_gtypes[ped_gtypes == -1] <- 2

  #remove the snps above
  ped_gtypes <- ped_gtypes[,-c(15,16)]
  expect_equal(ped_gtypes,show_genotypes(pop_b_gt_snp))

})

test_that("check summary stats are the same for gen_tibbles read in different ways",{
  geno_path <- system.file("extdata/pop_a.geno", package = "tidypopgen")
  pop_a_gt_geno <- gen_tibble(geno_path, quiet=TRUE, backingfile = tempfile(), valid_alleles = c("A","G","C","T"))

  bed_path <- system.file("extdata/pop_a.bed", package = "tidypopgen")
  pop_a_gt_bed <- gen_tibble(bed_path, quiet=TRUE, backingfile = tempfile(), valid_alleles = c("A","G","C","T"))

  ped_path <- system.file("extdata/pop_a.ped", package = "tidypopgen")
  pop_a_gt_ped <- gen_tibble(ped_path, quiet=TRUE, backingfile = tempfile(), valid_alleles = c("A","G","C","T"))

  # MAF
  bed_maf <- pop_a_gt_bed %>% loci_maf()
  ped_maf <- pop_a_gt_ped %>% loci_maf()
  geno_maf <- pop_a_gt_geno %>% loci_maf()
  expect_equal(bed_maf,ped_maf)
  expect_equal(bed_maf,geno_maf)
  expect_equal(ped_maf,geno_maf)

  # Missing loci
  bed_miss <- pop_a_gt_bed %>% loci_missingness()
  ped_miss <- pop_a_gt_ped %>% loci_missingness()
  geno_miss <- pop_a_gt_geno %>% loci_missingness()

  expect_equal(bed_miss,ped_miss)
  expect_equal(bed_miss,geno_miss)
  expect_equal(ped_miss,geno_miss)

  # Missing individuals
  bed_miss <- pop_a_gt_bed %>% indiv_missingness()
  ped_miss <- pop_a_gt_ped %>% indiv_missingness()
  geno_miss <- pop_a_gt_geno %>% indiv_missingness()

  expect_equal(bed_miss,ped_miss)
  expect_equal(bed_miss,geno_miss)
  expect_equal(ped_miss,geno_miss)

})


