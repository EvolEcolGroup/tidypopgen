test_indiv_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
test_genotypes <- rbind(c(2,2,0,1,1,2),
                        c(2,2,0,0,0,1),
                        c(2,2,0,0,1,1))
test_loci <- data.frame(name=paste0("rs",1:6),
                        chromosome=c(1,1,1,1,2,2),
                        position=c(3,5,65,343,23,456),
                        genetic_dist = as.integer(rep(0,6)),
                        allele_ref = c("A","T","C","G","C","T"),
                        allele_alt = c("T","C", NA,"C","G","A"))

test_gt <- gen_tibble(x = test_genotypes, indiv_meta = test_indiv_meta, loci = test_loci,
                      backingfile = tempfile())

test_that("ld clumping runs",{

  keep <- loci_ld_clump(test_gt, thr_r2 = 0.2, return_id=TRUE)
  expect_true(all.equal(keep, c(1, 2, 3, 4, 6)) == TRUE)

})



