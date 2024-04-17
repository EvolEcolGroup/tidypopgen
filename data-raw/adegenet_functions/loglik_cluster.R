# For the likelihood of cluster membership, we create
#
# binary matrices for 2 1 or 0 genotype
# AA
# AB
# BB
#
# matrices of frequences of allele A
# (frA) n loci x m clusters
#
# multiply
#
# 2 *AA * log(frA) + 2*B*B *log(1-frA) + AB *log(frA) + AB *log(1-frA) + AB * log(2)

## create a dataset for snapclust

# rows are loci, columns are clusters
freq<-rbind(c(0.3,0.05,0.4),
            c(0.1,0.8,0.1))


genotypes <- rbind(c(2,2),
                   c(1,0),
                   c(0,0),
                   c(1,1))

AA <-as.integer(genotypes==0)
AB <- as.integer(genotypes==1)
BB <- as.integer(genotypes==2)
freq_log2 <- matrix(log(2), ncol=ncol(freq), nrow=nrow(freq))
2 *(genotypes==0) %*% log(freq)+2*(genotypes==2) %*% log(1-freq)+
  (genotypes==1)%*%log(freq)+(genotypes==1)%*%log(1-freq)+(genotypes==1)%*%freq_log2

