

testthat::test_that(".genotypes_dot_prod computes correctly",{
  skip_if(TRUE)
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
   all.equal(res1,res2) # must be TRUE

   # genlight version
   library(adegenet)
   rownames(M) <- paste("ind", 1:100)
   x <- new("genlight",M)
   res11 <- glDotProd2(x,alleleAsUnit=FALSE, center=FALSE, scale=FALSE)


   for (i in 1:length(res1)){
     cat(i,"  ", identical(res1[[i]],res11[[i]]),"\n")
   }

   # not centred, not scaled
   res1 <- glDotProd(x,alleleAsUnit=FALSE, center=FALSE, scale=FALSE)
   res2 <- M %*% t(M)
   all.equal(res1,res2) # must be TRUE

   #  centred, not scaled
   res1 <- glDotProd(x,alleleAsUnit=FALSE, center=TRUE, scale=FALSE)
   M <- scalewt(M,center=TRUE,scale=FALSE)
   res2 <- M %*% t(M)
   all.equal(res1,res2) # must be TRUE

   #  centred, scaled
   res1 <- glDotProd(x,alleleAsUnit=FALSE, center=TRUE, scale=TRUE)
   M <- scalewt(M,center=TRUE,scale=TRUE)
   res2 <- M %*% t(M)
   all.equal(res1,res2) # must be TRUE

})
