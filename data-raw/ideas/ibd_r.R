# now test with missing data
test_na_gt <- gen_tibble(system.file("extdata/related/families.bed", package="tidypopgen"), quiet = TRUE,
                         backingfile = tempfile())


###########################

X_mat <- show_genotypes(test_na_gt)

## start by computing the expect probability of ibs given ibd
compute_EPrIBS_IBD <- function (X, CorrectFactor = TRUE){
  # matrix to store the expect probabilities
  EPrIBS_IBD <- matrix(0,ncol=3,nrow=3)
  nValid <- 0 # counter of valid snps for which we have estimates

  # iterate over each locus
  for (i in 1:ncol(X)){
    p <- sum(X[,i],na.rm = TRUE)/(sum(!is.na(X[,i]))*2)
    # Second, the expected probability of IBS i, given by IBD
    q = 1-p
    this_geno <- c(sum(X[,i]==0,na.rm = TRUE), sum(X[,i]==1, na.rm=TRUE),sum(X[,i]==2, na.rm=TRUE))
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
    if (!is.nan(a00) && !is.nan(a01) &&
        !is.nan(a02) && !is.nan(a11) && !is.nan(a12))
    {
      EPrIBS_IBD[1,1] = EPrIBS_IBD[1,1] + a00;
      EPrIBS_IBD[1,2] = EPrIBS_IBD[1,2] + a01;
      EPrIBS_IBD[1,3] = EPrIBS_IBD[1,3] + a02;
      EPrIBS_IBD[2,2] = EPrIBS_IBD[2,2] + a11;
      EPrIBS_IBD[2,3] = EPrIBS_IBD[2,3] + a12;
      nValid = nValid + 1
    }
  }

  EPrIBS_IBD[1,1] = EPrIBS_IBD[1,1]/ nValid;
  EPrIBS_IBD[2,1] = 0;
  EPrIBS_IBD[3,1] = 0;
  EPrIBS_IBD[1,2] = EPrIBS_IBD[1,2]/nValid;
  EPrIBS_IBD[2,2] = EPrIBS_IBD[2,2]/nValid;
  EPrIBS_IBD[3,2] = 0;
  EPrIBS_IBD[1,3] = EPrIBS_IBD[1,3]/ nValid;
  EPrIBS_IBD[2,3] = EPrIBS_IBD[2,3]/ nValid;
  EPrIBS_IBD[3,3] = 1
  return(EPrIBS_IBD)
}

EPrIBS_IBD <- compute_EPrIBS_IBD(X_mat)

## compute the ibs counts
X_mat0 <- X_mat==0
X_mat0[is.na(X_mat0)]<-0
X_mat1 <-X_mat==1
X_mat1[is.na(X_mat1)]<-0
X_mat2 <-X_mat==2
X_mat2[is.na(X_mat2)]<-0


## Should the snps be polarised to the major allele?
## otherwise it makes little sense to compare probabilities of 0, 1 and 2
ibs0 <- (X_mat0) %*% t(X_mat0)
ibs1 <- (X_mat1) %*% t(X_mat1)
ibs2 <- (X_mat2) %*% t(X_mat2)

n_ind <- nrow(X_mat)
# function to compute the quantities
Est_PLINK_Kinship <- function(IBS0, IBS1, IBS2, EPrIBS_IBD, constraint = FALSE) {
  # AM expected counts (p of each genotype times the total number of overlapping alleles)
  nIBS012 <- IBS0 + IBS1 + IBS2
  e00 <- EPrIBS_IBD[1, 1] * nIBS012
  e01 <- EPrIBS_IBD[1, 2] * nIBS012
  e11 <- EPrIBS_IBD[2, 2] * nIBS012
  e02 <- EPrIBS_IBD[1, 3] * nIBS012
  e12 <- EPrIBS_IBD[2, 3] * nIBS012
  e22 <- EPrIBS_IBD[3, 3] * nIBS012

  #browser()
  k0 <- IBS0 / e00
  k1 <- (IBS1 - k0 * e01) / e11
  k2 <- (IBS2 - k0 * e02 - k1 * e12) / e22

  # Bound IBD estimates to sum to 1, and fall within 0-1 range
  if (k0 > 1) { k0 <- 1; k1 <- k2 <- 0 }
  if (k1 > 1) { k1 <- 1; k0 <- k2 <- 0 }
  if (k2 > 1) { k2 <- 1; k0 <- k1 <- 0 }
  if (k0 < 0) { S <- k1 + k2; k1 <- k1 / S; k2 <- k2 / S; k0 <- 0 }
  if (k1 < 0) { S <- k0 + k2; k0 <- k0 / S; k2 <- k2 / S; k1 <- 0 }
  if (k2 < 0) { S <- k0 + k1; k0 <- k0 / S; k1 <- k1 / S; k2 <- 0 }

  if (constraint) {
    # Possibly constrain IBD estimates to within possible triangle
    # i.e. 0.5 0.0 0.5 is invalid
    #
    # Constraint : z1^2 - 4 z0 z2 >= 0
    #            : x^2 - 2 pi x + z2  = 0
    #              where pi = (z1 + 2 z2) / 2
    #
    # So the constaint can also be written as
    #              pi^2 >=  z2
    k2 <- 1 - k0 - k1
    pihat <- k1 / 2 + k2
    if (pihat^2 < k2) {
      k0 <- (1 - pihat) * (1 - pihat)
      k1 <- 2 * pihat * (1 - pihat)
    }
  }

  return(c(k0, k1, k2))
}

Est_PLINK_Kinship_pair <- function(x, ibs0, ibs1,ibs2, EPrIBS_IBD){
  Est_PLINK_Kinship(ibs0[x[1],x[2]],ibs1[x[1],x[2]],ibs2[x[1],x[2]], EPrIBS_IBD)
}

k0_mat <- matrix(0,ncol=n_ind, nrow=n_ind)
k1_mat <- matrix(0,ncol=n_ind, nrow=n_ind)
k2_mat <- matrix(0,ncol=n_ind, nrow=n_ind)

ibd_df <- Est_PLINK_Kinship(ibs0[1,2],ibs1[1,2],ibs2[1,2], EPrIBS_IBD)

indiv_comb <- combn(1:n_ind, 2)
foo <- t(apply(indiv_comb,2,Est_PLINK_Kinship_pair , ibs0, ibs1,ibs2, EPrIBS_IBD))
