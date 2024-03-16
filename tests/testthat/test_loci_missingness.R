test_that("snpbin_list_means computes correctly",{
  test_ind_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_genotypes <- rbind(c(1,1,0,1,1,2),
                          c(2,1,0,NA,0,NA),
                          c(2,2,0,0,1,NA))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=c(1,1,1,1,2,2),
                          position=c(3,5,65,343,23,456),
                          genetic_dist = as.integer(rep(0,6)),
                          allele_ref = c("a","t","c","g","c","t"),
                          allele_alt = c("t","c", NA,"c","g","a"))
  bed_path <- gt_write_bed_from_dfs(genotypes = test_genotypes,
                                    loci = test_loci,
                                    ind_meta = test_ind_meta,
                                    path_out = tempfile('test_data_'))
  test_gt <- gen_tibble(bed_path, quiet = TRUE)

  # na counts
  count_na <- function(x){sum(is.na(x))}
  n_na <- apply(test_genotypes, 2, count_na)
  expect_true(all(loci_missingness(test_gt$genotypes, as_counts = TRUE)==n_na))
  # convert to frequencies
  expect_true(all(loci_missingness(test_gt$genotypes)==n_na/nrow(test_genotypes)))
})
