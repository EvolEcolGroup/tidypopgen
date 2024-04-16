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
bed_path <- gt_write_bed_from_dfs(genotypes = test_genotypes,
                                  loci = test_loci,
                                  indiv_meta = test_indiv_meta,
                                  path_out = tempfile('test_data_'))
test_gt <- gen_tibble(bed_path, quiet = TRUE)


test_that("snp_ibs and gt_ibs computes ibs correctly",{
  test_fbm <- tidypopgen:::gt_get_bigsnp(test_gt)$genotypes
  test_ibs <- snp_ibs(test_fbm, as.counts=TRUE)
  # compare indiv 1 vs 2
  in_common<-sum(c(1,2,2,1,1,2))
  expect_identical(in_common, test_ibs$ibs[1,2])
  # check that we get the same result if we split the operation into two blocks
  test_ibs_2blocks <- snp_ibs(test_fbm, block.size = 3)
  expect_identical(test_ibs_2blocks$ibs[], test_ibs$ibs[])



  # now estimate it with gen_tibble
  test_ibs_gt <- gt_ibs(test_gt, as_counts = TRUE)
  expect_true(all.equal(test_ibs$ibs[], test_ibs_gt$ibs[],
                        check.attributes=FALSE))
})

