test_that("manually resorted table is respected by loci functions",{
  test_indiv_meta <- data.frame (id=c("a","b","c"),
                                 population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,0),
                          c(2,1,0,0,0,0),
                          c(2,2,0,0,1,1))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=paste0("chr",c(1,1,1,1,2,2)),
                          position=as.integer(c(3,5,65,343,23,456)),
                          genetic_dist = as.double(rep(0,6)),
                          allele_ref = c("A","T","C","G","C","T"),
                          allele_alt = c("T","C", NA,"C","G","A"))

  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)

  # test that loci functions can handle rearranged loci table
  original_freq_alt <- test_gt %>% loci_alt_freq()
  # create new tibble to reorder
  reorder_test_gt <- test_gt
  set.seed(123)
  new_order <- sample(1:nrow(show_loci(test_gt)))
  show_loci(reorder_test_gt) <- show_loci(test_gt)[new_order,]
  # estimate new frequencies
  reordered_freq_alt <- reorder_test_gt %>% loci_alt_freq()
  # check that the frequencies are the same
  expect_true(all(original_freq_alt[new_order] == reordered_freq_alt))
  # but clumping does generate an error
  expect_no_error(test_gt %>% loci_ld_clump(thr_r2 = 0.2))
  expect_error(reorder_test_gt %>% loci_ld_clump(thr_r2 = 0.2),
               "Your loci have been resorted")
})

test_that("gt_udpate_backingfile correctly updates",{
  test_indiv_meta <- data.frame (id=c("a","b","c","d","e"),
                                 population = c("pop1","pop1","pop2","pop2","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,0),
                          c(2,1,0,0,0,0),
                          c(2,2,0,0,1,1),
                          c(NA,2,0,NA,NA,0),
                          c(0,1,0,0,0,1))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=paste0("chr",c(2,1,1,1,1,2)),
                          position=as.integer(c(23,3,5,65,343,46)),
                          genetic_dist = as.double(rep(0,6)),
                          allele_ref = c("A","T","C","G","C","T"),
                          allele_alt = c("T","C", NA,"C","G","A"))
  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)
  # subset and reorder
  expect_false(is_loci_table_ordered(test_gt, error_on_false = FALSE))
  subset_reorder_test_gt <- test_gt %>% select_loci(c(2,4,5,1,6))
  expect_true(is_loci_table_ordered(subset_reorder_test_gt, error_on_false = FALSE))
  subset_reorder_test_gt <- subset_reorder_test_gt[c(2,1,4,5),]
  # now save the udpated backing matrix
  new_gt <- gt_udpate_backingfile(subset_reorder_test_gt, quiet = TRUE)
  # the new gt should be identical to the original one, minus the big indices
  expect_identical(show_genotypes(new_gt), show_genotypes(subset_reorder_test_gt))
  expect_identical(show_loci(new_gt)[,-1], show_loci(subset_reorder_test_gt)[,-1])
  # loci big index should now be sequential
  expect_identical(show_loci(new_gt)$big_index,1:count_loci(new_gt))
})

test_that("gt_order_loci reorders and regenerates backingfiles",{

  test_indiv_meta <- data.frame (id=c("a","b","c"),
                                 population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,0),
                          c(2,1,0,0,0,0),
                          c(2,2,0,0,1,1))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=paste0("chr",c(1,1,1,1,2,2)),
                          position=as.integer(c(3,5,65,343,23,456)),
                          genetic_dist = as.double(rep(0,6)),
                          allele_ref = c("A","T","C","G","C","T"),
                          allele_alt = c("T","C", NA,"C","G","A"))

  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)

  # create new tibble to reorder
  reorder_test_gt <- test_gt
  set.seed(123)
  new_order <- sample(1:nrow(show_loci(test_gt)))
  show_loci(reorder_test_gt) <- show_loci(test_gt)[new_order,]

  # write out to plink
  path <- paste0(tempfile(),"_outoforder")
  gt_as_plink(reorder_test_gt, path)
  path_bed <- paste0(path,".bed")

  # read back in
  reorder_test_gt <- gen_tibble(path_bed, quiet = TRUE)
  expect_true(is.integer(show_loci(reorder_test_gt)$chr_int))
  expect_true(is.integer(show_loci(reorder_test_gt)$position))

  # loci are out of order
  expect_false(is_loci_table_ordered(reorder_test_gt))

  # test gt_order_loci
  reorder_test_gt <- gt_order_loci(reorder_test_gt, use_current_table = FALSE, quiet = TRUE)

  # loci are now ordered, and backingfiles updated (v2)
  expect_true(is_loci_table_ordered(reorder_test_gt))
  expect_equal(normalizePath(gt_get_file_names(reorder_test_gt)),
               normalizePath(c(paste0(path,"_v2.rds"),paste0(path,"_v2.bk"))))

  # now save the gen_tibble
  gt_save(reorder_test_gt, quiet = TRUE)

  path_gt <- paste0(path,"_v2.gt")
  # reload and check the backingfile is updated
  reorder_test_gt_reload <- gt_load(path_gt)
  gt_get_file_names(reorder_test_gt_reload)

})

test_that("gt_order_loci use_current_table = TRUE",{

  test_indiv_meta <- data.frame (id=c("a","b","c"),
                                 population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,0),
                          c(2,1,0,0,0,0),
                          c(2,2,0,0,1,1))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=paste0("chr",c(2,1,1,1,1,2)),
                          position=as.integer(c(23,3,5,65,343,46)),
                          genetic_dist = as.double(rep(0,6)),
                          allele_ref = c("A","T","C","G","C","T"),
                          allele_alt = c("T","C", NA,"C","G","A"))

  path <- tempfile()
  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE, backingfile = path)

  # manually reorder the loci
  loci <- show_loci(test_gt)
  new_order <- c(2,3,4,5,1,6)
  show_loci(test_gt) <- show_loci(test_gt)[new_order,]
  expect_true(is_loci_table_ordered(test_gt))

  # test gt_order_loci use_current_table = TRUE
  test_gt <- gt_order_loci(test_gt, use_current_table = TRUE, quiet = TRUE)
  # now save the gen_tibble
  gt_save(test_gt, quiet = TRUE)

  # loci are still as they were reordered manually
  expect_equal(show_loci(test_gt)%>%select(-big_index),loci[new_order,]%>%select(-big_index))

  # but the backingfile has been updated
  expect_equal(normalizePath(gt_get_file_names(test_gt)),
               normalizePath(c(paste0(path,"_v2.rds"),paste0(path,"_v2.bk"))))

  # if use_current_table = TRUE, but the loci are not ordered, check for errors

  # manually reorder the loci for chromosomes to be out of order
  new_order <- c(5,1,2,3,4,6)
  show_loci(test_gt) <- show_loci(test_gt)[new_order,]

  # test gt_order_loci use_current_table = TRUE
  expect_error(gt_order_loci(test_gt, use_current_table = TRUE, quiet = TRUE),"All SNPs in a chromosome should be adjacent in the loci table")
  gt_save(test_gt, quiet = TRUE)

  # manually reorder the loci for positions to be out of order
  new_order <- c(3,4,2,5,1,6)
  show_loci(test_gt) <- show_loci(test_gt)[new_order,]
  expect_false(is_loci_table_ordered(test_gt))

  # test gt_order_loci use_current_table = TRUE
  expect_error(gt_order_loci(test_gt, use_current_table = TRUE, quiet = TRUE),"Your loci are not sorted within chromosomes")
})

