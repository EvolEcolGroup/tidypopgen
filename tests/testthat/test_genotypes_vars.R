testthat::test_that(".genotypes_vars computes correctly",{
  test_ind_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
  test_loci <- data.frame(name=paste0("rs",1:6),
                          chromosome=c(1,1,1,1,2,2),
                          position=c(3,5,65,343,23,432),
                          allele_ref = c("a","t","g","c","t","a"),
                          allele_alt = c("t","c","c","g","a","g"))

  # tests adapted from adegenet
  # same ploidy everywhere
  test_genotypes <- rbind(c(0,0,1,1,0, NA),
                          c(1,1,1,0,0,1),
                          c(0,0,0,1,1,1))
  test_gen <- gen_tibble(test_ind_meta, test_genotypes, test_loci,ploidy=c(1,1,1))
  f1 <- function(e) {return(mean((e-mean(e, na.rm=TRUE))^2, na.rm=TRUE))}
  expect_true(all.equal(.genotypes_vars(test_gen$genotypes),
                         apply(show_genotypes(test_gen, genotypes), 2, f1 )))
  expect_true(all.equal(.genotypes_vars(test_gen$genotypes,FALSE),
                         apply(show_genotypes(test_gen, genotypes), 2, f1 )))

  #  differences in ploidy
  test_genotypes <- rbind(c(0,0,1,1,0,NA), c(1,1,1,0,0,1), c(2,1,1,1,1,NA))
  test_gen <- gen_tibble(test_ind_meta, test_genotypes, test_loci,ploidy=c(1,1,2))
  temp <- sweep(show_genotypes(test_gen, genotypes), 1, c(1,1,2), "/")
  f2 <- function(e,w) {
    mu <- weighted.mean(e, w, na.rm=TRUE)
    res <- weighted.mean((e-mu)^2, w, na.rm=TRUE)
    return(res)
  }
  expect_true(all.equal(.genotypes_means(test_gen$genotypes),
                         apply(temp,2,weighted.mean, w=c(1,1,2), na.rm=TRUE)))
  expect_true(all.equal(.genotypes_vars(test_gen$genotypes), apply(temp, 2, f2,w=c(1,1,2) )))

  expect_true(all.equal(.genotypes_means(test_gen$genotypes,FALSE),
                         apply(temp,2,mean,na.rm=TRUE)))
  expect_true(all.equal(.genotypes_vars(test_gen$genotypes,FALSE), apply(temp,2,f1)))

})
