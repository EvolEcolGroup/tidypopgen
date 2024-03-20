# from https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS5.html

offspring.geno <- function(n.families, n.snps, fs = rep(0.5,n.snps), n.shared.parents=2){
  #INPUT:
  # n.families, number of families where each family produces two offspring (>0)
  # n.snps, number of independent SNPs used in simulation (>0)
  # fs, vector of allele 1 freqs for SNPs, length == n.snps, values >0 & <1
  # n.shared.parents, 0,1,2 shared parents for the two offspring in each family
  #OUTPUT:
  # X, the genotypes of (2*n.families) offspring, (2*n.families) x n.snps matrix with 0,1,2 entries

  stopifnot(n.families > 0)
  stopifnot(n.snps > 0)
  stopifnot(all(fs > 0 & fs < 1) & length(fs) == n.snps)
  stopifnot(n.shared.parents %in% 0:2)

  if(n.shared.parents == 2) parents = list(c(1,2), c(1,2)) #parents[[1]] are the parents of offspring 1
  if(n.shared.parents == 1) parents = list(c(1,2), c(3,2))
  if(n.shared.parents == 0) parents = list(c(1,2), c(3,4))
  n.parents = 4 - n.shared.parents #4, 3 or 2 for values of n.shared.parent of 0, 1 or 2

  X = matrix(0, nrow = 2*n.families, ncol = n.snps)
  for(ii in 1:n.families){ #each "family" means a pair of offspring that share 'n.shared.parents'
    x.parents = t(replicate(2*n.parents, rbinom(n.snps, size = 1, prob = fs) )) #2*n.parents parental genomes
    for(offs in 1:2){ #for two offspring within "family"
      #phase is the indicator of whether offs inherits each parents' 1st allele or not
      phase = t(replicate(2, rbinom(n.snps, size = 1, prob = 0.5))) #phase has one row for each parent
      for(i.parent in 1:2){
        for(ph in 0:1){
          loci = (phase[i.parent,] == ph) #which loci from i.parent have phase ph?
          #add to current offs' genotype i.parent's allele from the correct phase
          X[2*(ii-1) + offs, loci] =
            X[2*(ii-1) + offs, loci] + x.parents[2*parents[[offs]][i.parent] - ph, loci]
        }
      }
    }
  }
  return (X)
}
set.seed(123)
p = 10000 #SNPs
fs = runif(p, 0.2, 0.5) #MAF at each SNP is Uniform(0.2, 0.5)

X = rbind(offspring.geno(n.families = 5, n.snps = p, fs = fs, n.shared.parents = 0),
          offspring.geno(n.families = 5, n.snps = p, fs = fs, n.shared.parents = 1),
          offspring.geno(n.families = 5, n.snps = p, fs = fs, n.shared.parents = 2))

X = X[,apply(X,2,var) > 0] #remove possible monomorphic variants

# IBS calculation
IBS.2 = ( (X==2) %*% t(X==2) + (X==1) %*% t(X==1) + (X==0) %*% t(X==0) )/p #prop. of loci with same genotypes
IBS.1 = ( (X==1) %*% t(X==0 | X==2) + (X==0 | X==2) %*% t(X==1) )/p #prop. of loci with one allele shared
IBS = IBS.2 + 0.5*IBS.1 #prop. of genome IBS


X_a <- X[,1:5000]
X_b <- X[,5001:10000]

# IBS2_fun <- function (X){
#   ( (X==2) %*% t(X==2) + (X==1) %*% t(X==1) + (X==0) %*% t(X==0) )
# }
#
# foo<- IBS2_fun(X_sub1)
# foo2 <- IBS2_fun(X_sub2)

#identical((foo+foo2)/p, IBS.2)


# break it down
same_2 <- function(X){
  (X==2) %*% t(X==2) + (X==1) %*% t(X==1) + (X==0) %*% t(X==0)
}


# IBS calculation
IBS.2.c = ( (X==2) %*% t(X==2) + (X==1) %*% t(X==1) + (X==0) %*% t(X==0) ) #prop. of loci with same genotypes
IBS.1.c = ( (X==1) %*% t(X==0 | X==2) + (X==0 | X==2) %*% t(X==1) ) #prop. of loci with one allele shared
IBS.c = 2*IBS.2.c + IBS.1.c #prop. of genome IBS

all.equal(IBS.c/(2*p),IBS)


ibs_counts <- function(X,ind){
  X_sub <- X[,ind] # get a slice of the matrix
  2 * ( (X_sub==2) %*% t(X_sub==2) + (X_sub==1) %*% t(X_sub==1) + (X_sub==0) %*% t(X_sub==0) ) + #prop. of loci with same genotypes
    ( (X_sub==1) %*% t(X_sub==0 | X_sub==2) + (X_sub==0 | X_sub==2) %*% t(X_sub==1) ) #prop. of loci with one allele shared
}

all.equal(IBS_counts(X,1:5000)+IBS_counts(X,5001:10000),IBS_counts(X,1:ncol(X)))

IBS_prop <-IBS_counts / n_loci


### Now let's try to implement with FBM
geno <- bigstatsr::FBM.code256(
  nrow = nrow(X),
  ncol = ncol(X),
  code = bigsnpr::CODE_012,
  backingfile = tempfile("geno"),
  init = NULL,
  create_bk = TRUE
)
geno[] <- X



my_counts <- bigstatsr::big_apply(geno, a.FUN = ibs_counts, a.combine = "plus")
all.equal(my_counts/(2*p),IBS)


snp_ibs_R <- function(X,
                    row.ind = bigstatsr::rows_along(X),
                    col.ind = bigstatsr::cols_along(X),
                    as.counts = TRUE,
                    block.size = block_size(nrow(X), ncores)) {
  # TODO check that the dataset has no NAs!!!!

  # IBS (as counts) for a matrix slice
  # TODO this is prime material for converting with RccpArmadillo
  ibs_counts <- function(X, ind){
    X_sub <- X[,ind] # get a slice of the matrix
    X_sub_t <- t(X_sub)
    2 * ( (X_sub==2) %*% (X_sub_t==2) + (X_sub==1) %*% (X_sub_t==1) + (X_sub==0) %*% (X_sub_t==0) ) + #prop. of loci with same genotypes
      ( (X_sub==1) %*% (X_sub_t==0 | X_sub_t==2) + (X_sub==0 | X_sub==2) %*% (X_sub_t==1) ) #prop. of loci with one allele shared
  }
  ibs_counts_matrix <- bigstatsr::big_apply(X,
                                            a.FUN = ibs_counts,
                                            ind = col.ind,
                                            a.combine = "plus",
                                            block.size = block.size)
  if (as_counts){
    return(ibs_counts_matrix)
  } else {
    ibs_counts_matrix/(2*length(col.ind))
  }
}

ibs_with_big <- snp_ibs(geno)
all.equal(ibs_with_big, IBS)





############################
# some fun with matrices
############################
# create a copy
geno2 <- geno$copy()

# the codes
# create code for X==0
code_0 <- rep(0, 256)
code_0[1:3] <- c(1,0,0)
# create code for X==1
code_1 <- rep(0, 256)
code_1[1:3] <- c(0,1,0)
# create code for X==2
code_2 <- rep(0, 256)
code_2[1:3] <- c(0,0,1)

geno <- add_code256(geno, code_1)
first_part <- big_tcrossprodSelf(geno)
geno <- add_code256(geno, code_0)
geno2 <- add_code256(geno, code_2)
second_part <- big_prodMat(geno,t(geno2[]))
third_part <- big_prodMat(geno2,t(geno[]))

identical(first_part[], (X==1) %*% t(X==1))
identical(second_part[], (X==0) %*% t(X==2))
identical(third_part[], (X==2) %*% t(X==0))


k <- 2*first_part[] -2*(second_part[]+third_part[])
k2 <- 2*((X==1) %*% t(X==1) - 2*((X==0) %*% t(X==2) + (X==2) %*% t(X==0)) )

denominator = matrix(rep(rowSums(X==1), nrow(X)), nrow = nrow(X), byrow = T) +
  matrix(rep(rowSums(X==1), nrow(X)), nrow = nrow(X), byrow = F)
king.r = 2*((X==1) %*% t(X==1) - 2*((X==0) %*% t(X==2) + (X==2) %*% t(X==0)) ) / denominator
#(suppressing commands of matrix printing from output)



