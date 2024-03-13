test_that("select_loci subsets correctly",{
  test_ind_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,2),
                          c(2,1,0,NA,0,NA),
                          c(2,2,0,0,1,NA))
  test_loci <- data.frame(name=c(paste0("rs",1:4),paste0("x",1:2)),
                          chromosome=c(1,1,1,1,2,2),
                          position=c(3,5,65,343,23,456),
                          genetic_dist = as.integer(rep(0,6)),
                          allele_ref = c("a","t","c","g","c","t"),
                          allele_alt = c("t","c", NA,"c","g","a"))
  file_plink <- tempfile('test_data_')
  make_test_bed(test_genotypes, test_loci, test_ind_meta, file_plink)
  file_plink<-paste0(file_plink,".bed")
  # convert bed to bigsnp
  path_rds <- bigsnpr::snp_readBed(file_plink, backingfile = tempfile("test_bigfile_"))
  # convert to gen_tibble
  test_gt <- gen_tibble(path_rds)

  # select snps with an rs
  test_gt_sub <- test_gt %>% select_loci (starts_with("rs"))
  expect_true(!any(c("x1","x2") %in% show_loci_names(test_gt_sub)))
  # subsetting by id with reordering
  test_gt_sub <- test_gt %>% select_loci (c(3,1,5))
  expect_identical(c("rs3","rs1","x1"), show_loci_names(test_gt_sub))
  # get everything
  test_gt_sub <- test_gt %>% select_loci (everything())
  expect_identical(test_gt, test_gt_sub)
  # use 2:4 range expressions
  test_gt_sub <- test_gt %>% select_loci (2:4)
  expect_identical(test_loci$name[2:4], show_loci_names(test_gt_sub))

})
