test_that("indiv_missingness computes correctly",{
  test_indiv_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,2),
                          c(2,1,0,NA,0,NA),
                          c(2,2,0,0,1,NA))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=c(1,1,1,1,2,2),
                          position=c(3,5,65,343,23,456),
                          genetic_dist = as.integer(rep(0,6)),
                          allele_ref = c("A","T","C","G","C","T"),
                          allele_alt = c("T","C", NA,"C","G","A"))

  bed_path <- gt_write_bed_from_dfs(genotypes = test_genotypes,
                                    loci = test_loci,
                                    indiv_meta = test_indiv_meta,
                                    path_out = tempfile('test_data_'))
  test_gt <- gen_tibble(bed_path, quiet = TRUE)

  sum_na <- function(x){sum(is.na(x))}
  # feeding the genotypes directly
  expect_true(all(indiv_missingness(test_gt$genotypes)==
                    apply(test_genotypes,1,sum_na)/ncol(test_genotypes)))
  # passing tibble
  expect_true(all(indiv_missingness(test_gt)==
                    apply(test_genotypes,1,sum_na)/ncol(test_genotypes)))
})
