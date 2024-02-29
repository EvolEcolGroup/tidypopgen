

testthat::test_that(".genotypes_dot_prod computes correctly",{
   # tests taken from adegenet
   # TESTING DOT PRODUCTS
   set.seed(123)
   M <- matrix(sample(c(0,1), 100*1e3, replace=TRUE), nrow=100)
   n_loci <- ncol(M)
   test_loci <- data.frame(name=paste0("rs",n_loci),
                           chromosome=rep(1,n_loci),
                           position=rep(0,n_loci),
                           allele_ref =rep("a",n_loci),
                           allele_alt = rep("t",n_loci))
   test_ind_meta <- data.frame (id=paste0("ind",nrow(M)),
                                population = rep(c("pop1","pop2"),each=50))

   test_gen <- gen_tibble(test_ind_meta, M, test_loci, ploidy = rep(1,100))


   # not centred, not scaled
   res1 <- .genotypes_dot_prod(test_gen$genotypes,alleles_as_units=FALSE, center=FALSE, scale=FALSE)
   res2 <- M %*% t(M)
   expect_true(all.equal(res1, res2, check.attributes=FALSE))

   #  centred, not scaled
   res1c <- .genotypes_dot_prod(test_gen$genotypes,alleles_as_units=FALSE, center=TRUE, scale=FALSE)
   M1 <- scalewt(M,center=TRUE,scale=FALSE)
   res2c <- M1 %*% t(M1)
   expect_true(all.equal(res1c,res2c, check.attributes=FALSE))

   #  centred, scaled
   res1cs <- .genotypes_dot_prod(test_gen$genotypes,alleles_as_units=FALSE, center=TRUE, scale=TRUE)
   M2 <- scalewt(M,center=TRUE,scale=TRUE)
   res2cs <- M2 %*% t(M2)
   expect_true(all.equal(res1cs,res2cs, check.attributes=FALSE))

   #TODO add tests for parallel version, testing against the serial objects


})
