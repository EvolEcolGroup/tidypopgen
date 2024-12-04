# create file
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

# create a gen_tibble with individuals and loci subsetted and reordered
# create file
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
expect_false(is_loci_table_ordered(test_gt))
subset_reorder_test_gt <- test_gt %>% select_loci(c(2,4,5,1,6))
expect_true(is_loci_table_ordered(subset_reorder_test_gt))
subset_reorder_test_gt <- subset_reorder_test_gt[c(2,1,4,5),]
# now save the udpated backing matrix
new_gt <- gt_update_backing_file(subset_reorder_test_gt)
# the new gt should be identical to the original one, minus the big indeces
expect_identical(show_genotypes(new_gt), show_genotypes(subset_reorder_test_gt))
expect_identical(show_loci(new_gt)[,-1], show_loci(subset_reorder_test_gt)[,-1])
# loci big index should now be sequential
expect_identical(show_loci(new_gt)$big_index,1:count_loci(new_gt))
