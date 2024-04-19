# ibd with matrix operations

# dataset
test_indiv_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
test_genotypes <- rbind(c(1,1,0,1,0,0),
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

###########################

X <- test_genotypes

ibs0 <- (X==0) %*% t(X==0)
ibs1 <- (X==1) %*% t(X==1)
ibs2 <- (X==2) %*% t(X==2)

genos <- table(X)
n_valid_alleles <- ibs0 + ibs1 + ibs2


## start by computing the expect probability of ibd

CorrectFactor <- TRUE
nValid <- 0 # counter of valid snps for which we have estimates

# iterate over each locus
for (i in 1:ncol(X)){
  p <- sum(X[,i])/(nrow(X)*2)
  # Second, the expected probability of IBS i, given by IBD
  q = 1-p
  this_geno <- c(sum(X[,i]==0), sum(X[,i]==1),sum(X[,i]==2))
  Na = n = sum(this_geno)
  x = 2*this_geno[1] + this_geno[2] ## AA AB
  y = 2*this_geno[3] + this_geno[2] ## BB AB
  a00 =a01 =a02=a11=a12 = 0

  if (CorrectFactor)  {
    a00 =
      2*p*p*q*q * ( (x-1)/x * (y-1)/y *
                      (Na/(Na-1)) *(Na/(Na-2)) * (Na/(Na-3)) );
    a01 =
      4*p*p*p*q * ( (x-1)/x * (x-2)/x * (Na/(Na-1)) *
                      (Na/(Na-2)) * (Na/(Na-3)) ) + 4*p*q*q*q * ( (y-1)/y * (y-2)/y *
                                                                    (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
    a02 =
      q*q*q*q * ( (y-1)/y * (y-2)/y * (y-3)/y *
                    (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) ) + p*p*p*p *
      ( (x-1)/x * (x-2)/x * (x-3)/x * (Na/(Na-1)) * (Na/(Na-2)) *
          (Na/(Na-3)) ) + 4*p*p*q*q * ( (x-1)/x * (y-1)/y * (Na/(Na-1)) *
                                          (Na/(Na-2)) * (Na/(Na-3)) );
    a11 =
      2*p*p*q * ( (x-1)/x *  Na/(Na-1) * Na/(Na-2) ) +
      2*p*q*q * ( (y-1)/y *  Na/(Na-1) * Na/(Na-2) );
    a12 =
      p*p*p * ((x-1)/x * (x-2)/x *  Na/(Na-1) * Na/(Na-2)) +
      q*q*q * ( (y-1)/y * (y-2)/y *  Na/(Na-1) * Na/(Na-2)) +
      p*p*q * ( (x-1)/x * Na/(Na-1) * Na/(Na-2) ) +
      p*q*q * ((y-1)/y  * Na/(Na-1) * Na/(Na-2));
  } else {
    a00 = 2*p*p*q*q;
    a01 = 4*p*p*p*q + 4*p*q*q*q;
    a02 = q*q*q*q + p*p*p*p + 4*p*p*q*q;
    a11 = 2*p*p*q + 2*p*q*q;
    a12 = p*p*p + q*q*q + p*p*q + p*q*q;
  }

  # now check that a00, a.. all finite, then add them
  # this is only an issue when we use the correction factor
  # alternatively, couldn't we check that Na>3 and x and y >1???
  # it should be much faster than check for finite values

  # TO TRANSLATE INTO R
  if (R_FINITE(a00) && R_FINITE(a01) &&
      R_FINITE(a02) && R_FINITE(a11) && R_FINITE(a12))
  {
    EPrIBS_IBD[0][0] += a00;
    EPrIBS_IBD[0][1] += a01;
    EPrIBS_IBD[0][2] += a02;
    EPrIBS_IBD[1][1] += a11;
    EPrIBS_IBD[1][2] += a12;
    nValid++;
  }
}

EPrIBS_IBD[0][0] /= nValid;
EPrIBS_IBD[1][0] = 0;
EPrIBS_IBD[2][0] = 0;
EPrIBS_IBD[0][1] /= nValid;
EPrIBS_IBD[1][1] /= nValid;
EPrIBS_IBD[2][1] = 0;
EPrIBS_IBD[0][2] /= nValid;
EPrIBS_IBD[1][2] /= nValid;
EPrIBS_IBD[2][2] = 1;


}
