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


# this also tests show_genotypes and show_loci
test_that("write a bed file",{
  bed_path <- gt_as_plink(test_gt, file = paste0(tempfile(),".bed"))
  # now read the file back in
  test_gt2 <- gen_tibble(bed_path, quiet=TRUE)
  ## continue here

  # because of the different backing file info, we cannot use identical on the whole object
  expect_true(identical(show_genotypes(test_gt), show_genotypes(test_gt2)))
  expect_true(identical(test_gt %>% select(-genotypes),
                        test_gt2 %>% select(-genotypes)))

  #check gt_as_plink converts the NA missing allele to 0
  expect_true(is.na(show_loci(test_gt2)$allele_alt[3]))


  # now write it as a ped
  ped_path <- gt_as_plink(test_gt, file = paste0(tempfile(),".ped"), type = "ped")
  test_gt3 <- gen_tibble(ped_path, quiet=TRUE)
  # the gen tibble from the bed and ped should contain the same information
  expect_true(all.equal(show_loci(test_gt3),show_loci(test_gt2), check.attributes=FALSE))
  expect_true(all.equal(show_genotypes(test_gt3),show_genotypes(test_gt2)))

  # write it as raw
  raw_path <- gt_as_plink(test_gt,file=tempfile(),type="raw")
  # check why sex has been lost (it was not read in properly?!?)
  # look at line 119 of gen_tibble
  #tools::md5sum(raw_path)==tools::md5sum(system.file("extdata","pop_a.raw",package="tidypopgen"))

})

