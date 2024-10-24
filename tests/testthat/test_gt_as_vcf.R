# unit test for gt_as_vcf
# create file
test_indiv_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
test_genotypes <- rbind(c(1,1,0,1,1,0),
                        c(2,1,0,0,0,0),
                        c(2,2,0,0,NA,1))
test_loci <- data.frame(name=paste0("rs",1:6),
                        chromosome=paste0("chr",c(1,1,1,1,2,2)),
                        position=as.integer(c(3,5,65,343,23,456)),
                        genetic_dist = as.integer(rep(0,6)),
                        allele_ref = c("A","T","C","G","C","T"),
                        allele_alt = c("T","C", NA,"C","G","A"))

test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)

test_that("test we can write and read a vcf",{
  vcf_path <- gt_as_vcf(test_gt, file = paste0(tempfile(),".vcf"))
  # now read the file back in
  test_gt2 <- gen_tibble(vcf_path, quiet=TRUE, parser = "vcfR", backingfile = tempfile())
  # because of the different backing file info, we cannot use identical on the whole object
  expect_true(identical(show_genotypes(test_gt), show_genotypes(test_gt2)))
  #check gt_as_plink converts the NA missing allele to 0
  expect_true(is.na(show_loci(test_gt2)$allele_alt[3]))
  # now use the cpp parser
  test_gt2 <- gen_tibble(vcf_path, quiet=TRUE, parser = "cpp", backingfile = tempfile())

  expect_true(identical(show_genotypes(test_gt), show_genotypes(test_gt2)))
  #check gt_as_plink converts the NA missing allele to 0
  expect_true(is.na(show_loci(test_gt2)$allele_alt[3]))
})

test_that("test reading and writing is equivalent",{
  # read in .vcf
  vcf_path <- system.file("extdata/pop_a.vcf", package = "tidypopgen")
  pop_a_vcf <- gen_tibble(vcf_path, quiet=TRUE,backingfile = tempfile(), parser="vcfR")
  # write it out again
  file <- paste0(tempfile(),".vcf")
  write_vcf_path <- gt_as_vcf(pop_a_vcf, file = file)
  # reread
  pop_a_vcf_rewrite <- gen_tibble(write_vcf_path, quiet=TRUE,backingfile = tempfile(), parser="vcfR")
  # because of the different backing file info, we cannot use identical on the whole object
  expect_true(identical(show_genotypes(pop_a_vcf), show_genotypes(pop_a_vcf_rewrite)))
  expect_true(identical(show_loci(pop_a_vcf), show_loci(pop_a_vcf_rewrite)))
})

#TESTS TO SHOW PROBLEM WITH CHUNKING:

# when read file back in with cpp it lets you read it in despite duplicates so can
# see the issue
test_that("test reading and writing with chunking is equivalent cpp",{
  # write out to vcf with chunk size - smaller chunk size = more duplicates
  vcf_path_chunked <- gt_as_vcf(test_gt, file = paste0(tempfile(),".vcf"), chunk_size = 2)
  # read vcf back in with cpp
  test_gt_chunked <- gen_tibble(vcf_path_chunked, quiet=TRUE, parser = "cpp", backingfile = tempfile("anolis_chunk_"))
  # check they are the same
  expect_true(identical(show_genotypes(test_gt), show_genotypes(test_gt_chunked)))
  expect_true(identical(show_loci(test_gt), show_loci(test_gt_chunked)))
  expect_equal(count_loci(test_gt_chunked), count_loci(test_gt))
})

# when read file back in with vcfR it throws an error because of the duplicates
# so wont read in
test_that("test reading and writing with chunking is equivalent rvcf",{
  # write out to vcf with chunk size
  vcf_path_chunked <- gt_as_vcf(test_gt, file = paste0(tempfile(),".vcf"), chunk_size = 2)
  # read vcf back in with vcfR
  test_gt_chunked <- gen_tibble(vcf_path_chunked, quiet=TRUE, parser = "vcfR", backingfile = tempfile("anolis_chunk_"))
  # check they are the same
  expect_true(identical(show_genotypes(test_gt), show_genotypes(test_gt_chunked)))
  expect_true(identical(show_loci(test_gt), show_loci(test_gt_chunked)))
  expect_equal(count_loci(test_gt_chunked), count_loci(test_gt))
})








