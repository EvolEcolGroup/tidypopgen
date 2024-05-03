library(hierfstat)
set.seed(123)
dos<-matrix(sample(0:2,size=10000,replace=TRUE),ncol=100)
pop=rep(1:5,each=20)
fs.dosage(dos,pop=pop)

dos_match <- matching(dos)

####
test_genotypes <- rbind(c(1,1,0,1,1,0),
                        c(2,1,0,NA,0,0),
                        c(2,NA,0,0,1,1))
test_indiv_meta <- data.frame (id=c("a","b","c"),
                               population = c("pop1","pop1","pop2"))
test_loci <- data.frame(name=paste0("rs",1:6),
                        chromosome=paste0("chr",c(1,1,1,1,2,2)),
                        position=as.integer(c(3,5,65,343,23,456)),
                        genetic_dist = as.integer(rep(0,6)),
                        allele_ref = c("A","T","C","G","C","T"),
                        allele_alt = c("T","C", NA,"C","G","A"))

test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)



dos_hier_match <- hierfstat::matching(test_genotypes)

test_fbm <- tidypopgen:::.gt_get_bigsnp(test_gt)$genotypes
test_as <- snp_allele_sharing(test_fbm)



dos <- test_genotypes

na <- matrix(rep(1,prod(dim(dos))),ncol=ncol(dos))
ina<-which(is.na(dos))
na[ina]<-0
dos[ina]<-1
Mij <- 1/2 * (1+1/tcrossprod(na) * tcrossprod(dos-1))
# Mij is the allele sharing matrix
