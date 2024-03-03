a1<-15
n1<-25
# expect het for a locus
h1 <- (a1*(n1-a1))/(n1*(n1-1))
h1



#########
a1 <- colSums2(as.matrix(pop1),na.rm=T)
a2 <- colSums2(as.matrix(pop2),na.rm=T)
n1 <- apply(as.matrix(pop1),2,function(x) 2*sum(!is.na(x)))
n2 <- apply(as.matrix(pop2),2,function(x) 2*sum(!is.na(x)))

h1 <- (a1*(n1-a1))/(n1*(n1-1))
h2 <- (a2*(n2-a2))/(n2*(n2-1))

N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
D <- N + h1 + h2

F <- sum(N, na.rm=T)/sum(D, na.rm=T)


pairwise_fst <- function(.x) {
  # check that the data are diploid

  #@TODO
  # loci sums
  # do something clever like first count na and then apply addition by column
}

# make such a function using loci_freq as a template
# then do the same for het_obs
# finally use these two functions for ind functions (by taking the mean)
loci_het_exp <- function(.x){
  sums <- snpbin_list_sums(.x, alleles_as_units = TRUE)
  n <- snpbin_list_n(.x, alleles_as_units = TRUE)
  return((sums*(n-sums))/(n*(n-1)))
}

pairwise_fst <- function(.x, .x2){

}
