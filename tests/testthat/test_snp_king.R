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


# this also tests show_genotypes and show_loci
test_that("snp_king_r and gt_king compute king-robust correctly",{
  test_fbm <- tidypopgen:::gt_get_bigsnp(test_gt)$genotypes
  test_king <- snp_king(test_fbm)
  # king by hand
  # code from https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS5.html
  X <- test_genotypes
  denominator = matrix(rep(rowSums(X==1), nrow(X)), nrow = nrow(X), byrow = T) +
    matrix(rep(rowSums(X==1), nrow(X)), nrow = nrow(X), byrow = F)
  king.r = 2*((X==1) %*% t(X==1) - 2*((X==0) %*% t(X==2) + (X==2) %*% t(X==0)) ) / denominator
  expect_identical(king.r, test_king)
  # check that we get the same result if we split the operation into two blocks
  test_king_2blocks <- snp_king(test_fbm, block.size = 3)
  expect_identical(test_king_2blocks, test_king)

  # now estimate it with gen_tibble
  test_king_gt <- gt_king(test_gt)
  expect_true(all.equal(test_king, test_king_gt,
                        check.attributes=FALSE))
})

