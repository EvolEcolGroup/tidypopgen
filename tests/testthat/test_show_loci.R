test_that("show_loci gets and sets information",{
  test_ind_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,2),
                          c(2,1,0,NA,0,NA),
                          c(2,2,0,0,1,NA))
  test_loci <- data.frame(name=c(paste0("rs",1:4),paste0("x",1:2)),
                          chromosome=as.integer(c(1,1,1,1,2,2)),
                          position=as.integer(c(3,5,65,343,23,456)),
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
  # check that we retrieve the info we put in (as a tibble)
  expect_identical(show_loci(test_gt) %>% select(-big_index),as_tibble(test_loci))
  # now change it directly on the genotype column
  test_loci2 <- test_loci %>% dplyr::mutate(chromosome = "new")
  show_loci(test_gt$genotypes) <- test_loci2
  expect_identical(show_loci(test_gt), as_tibble(test_loci2))
  test_loci3 <- test_loci %>% dplyr::mutate(chromosome = "newer")
  show_loci(test_gt) <- test_loci3
  expect_identical(show_loci(test_gt), as_tibble(test_loci3))
  # with some proper dplyr
  show_loci(test_gt) <- show_loci(test_gt) %>% mutate(chromosome="old")
  expect_true(all(show_loci(test_gt)$chromosome=="old"))
  test_loci3 <- test_loci3[-1,]
  expect_error(show_loci(test_gt) <- test_loci3,
               "the replacement loci tibble does")
})
