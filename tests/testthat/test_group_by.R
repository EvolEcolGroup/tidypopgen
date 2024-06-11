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
test_that("inheritance of group tibbles",{
  # test operations on standard gen_tibbles
  expect_true(inherits(test_gt, "gen_tbl"))
  # add column
  test_gt <- test_gt %>% mutate(x=c(1,2,3))
  expect_true(inherits(test_gt, "gen_tbl"))
  # remove column
  test_sub_gt <- test_gt %>% select(-x)
  expect_true(inherits(test_sub_gt, "gen_tbl"))
  # remove genotype column
  test_sub_gt <- test_gt %>% select(-genotypes)
  expect_false(inherits(test_sub_gt, "gen_tbl"))
  # filter
  test_sub_gt <- test_gt %>% filter(population=="pop1")
  expect_true(inherits(test_sub_gt, "gen_tbl"))

  # now group it
  test_group_gt <- test_gt %>% group_by(population)
  expect_true(inherits(test_group_gt, "gen_tbl"))
  expect_true(inherits(test_group_gt, "grouped_gen_tbl"))
  # slice rows
  test_group_gt <- test_group_gt %>% filter(population == "pop1")
  expect_true(inherits(test_group_gt, "gen_tbl"))
  expect_true(inherits(test_group_gt, "grouped_gen_tbl"))
  # add columns
  test_group_gt <- test_group_gt %>% mutate(n=c(1,2))
  expect_true(inherits(test_group_gt, "gen_tbl"))
  expect_true(inherits(test_group_gt, "grouped_gen_tbl"))

  test_group <- test_group_gt %>% select(-genotypes)
  expect_false(inherits(test_group, "gen_tbl"))
  expect_false(inherits(test_group, "grouped_gen_tbl"))

  test_ungroup_gt <- ungroup(test_group_gt)
  expect_true(inherits(test_ungroup_gt, "gen_tbl"))
})
