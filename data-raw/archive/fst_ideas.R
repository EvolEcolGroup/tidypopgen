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

pop_pairwise_fst <- function(.x, by_locus=FALSE){

  warning("this function is not properly tested yet!!!")
  # check matrix(unlist(z, use.names = FALSE), ncol = 10, byrow = TRUE)
  # is known to be faster than do.call(rbind,f)
  # see https://stackoverflow.com/questions/13224553/how-to-convert-a-huge-list-of-vector-to-a-matrix-more-efficiently
  # do not modify the approach used below:

  n_loci <- adegenet::nLoc(.x$genotypes[[1]])
  # sum alt alleles over each locus in each group
  sums <- matrix(unlist(.x %>%
                          group_map(.f=~snpbin_list_sums(.x$genotypes, alleles_as_units = TRUE)),
                        use.names = FALSE), ncol = n_loci, byrow = TRUE)
  # get the total  number of alleles (i.e. removing NAs) for each locus in each group
  n <- matrix(unlist(.x %>%
                       group_map(.f=~snpbin_list_n(.x$genotypes, alleles_as_units = TRUE)),
                        use.names = FALSE), ncol = n_loci, byrow = TRUE)
  # function to compute het by row
  het_exp_by_row <- function(i, sums, n){(sums[i,]*(n[i,]-sums[i,]))/(n[i,]*(n[i,]-1))}
  # get het at each locus for each population
  het <- matrix(unlist(lapply(1:nrow(sums), het_exp_by_row, sums, n),
                     use.names = FALSE), ncol = n_loci, byrow = TRUE)

  # get the grouping column, and creat all pairwise combination of indices
  .group_levels = .x %>% group_keys()
  pairwise_combn <- t(combn(nrow(.group_levels),2))
  numerator <- matrix(NA_real_, nrow = nrow(pairwise_combn), ncol = n_loci)
  denominator <- matrix(NA_real_, nrow = nrow(pairwise_combn), ncol = n_loci)
  for (i_row in seq_len(nrow(pairwise_combn))){
    pop1 <- pairwise_combn[i_row,1]
    pop2 <- pairwise_combn[i_row,2]
    numerator[i_row,] <- (sums[pop1,]/n[pop1,] - sums[pop2,]/n[pop2,])^2 -
      het[pop1,]/n[pop1,] - het[pop2,]/n[pop2,]
    denominator[i_row,] <- numerator[i_row,] + het[pop1,] + het[pop2,]
  }

  # if we want whole genome level estimates
  fst <- tibble(!!paste0(names(.group_levels),"_",1) := .group_levels[pairwise_combn[,1],1] %>% pull(1),
                !!(paste0(names(.group_levels),"_",2)) := .group_levels[pairwise_combn[,2],1] %>% pull(1))
  if (!by_locus){
    fst_wg_by_row <- function(i, numerator, denominator){
      return(sum(numerator[i,],na.rm=TRUE)/sum(denominator[i,],na.rm=TRUE))
    }
    fst <- fst %>% mutate(value = unlist(lapply(seq_len(nrow(fst)),fst_wg_by_row, numerator, denominator)))

    return(fst)
  } else {
    # if we want estimates by locus
    fst_by_row <- function(i, numerator, denominator){
      return(numerator[i,]/denominator[i,])
    }
    fst_vals <- matrix(unlist(lapply(seq_len(nrow(numerator)),fst_by_row, numerator, denominator),
                         use.names = FALSE), ncol = n_loci, byrow = TRUE)
    colnames(fst_vals) <- show_loci_names(.x)
    return(fst %>% cbind(fst_vals))
  }
}
