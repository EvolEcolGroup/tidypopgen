test_genotypes <- rbind(c(1,1,0,1,1,0),
                        c(2,1,0,NA,0,0),
                        c(2,NA,0,0,1,1),
                        c(1,0,0,1,0,0),
                        c(1,2,0,1,2,1),
                        c(0,0,0,0,NA,1),
                        c(0,1,1,0,1,NA))
test_indiv_meta <- data.frame (id=c("a","b","c","d","e","f","g"),
                               population = c("pop1","pop1","pop2","pop2","pop1","pop3","pop3"))
test_loci <- data.frame(name=paste0("rs",1:6),
                        chromosome=paste0("chr",c(1,1,1,1,2,2)),
                        position=as.integer(c(3,5,65,343,23,456)),
                        genetic_dist = as.double(rep(0,6)),
                        allele_ref = c("A","T","C","G","C","T"),
                        allele_alt = c("T","C", NA,"C","G","A"))

test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)
test_gt <- test_gt %>% group_by(population)

# a developer function to create various count summaries of a population, used to
# compute more complex statistics (e.g. pairwise fst, etc.).
# Specifically, we compute:
# sums_alt sum alleles for alt
# sums_ref sums alleles for ref
# n sums of all alleles
# het_obs observed heterozygosity
# het_exp expected heterozygosity
.gt_pop_freqs <- function(.x){
  counts <- bigstatsr::big_counts( .gt_get_bigsnp(.x)$genotypes,
                                   ind.row =.gt_bigsnp_rows(.x),
                                   ind.col = .gt_bigsnp_cols(.x))
  sums_alt <- apply(counts,2,function(x) x[2]+2*x[3])
  n <- apply(counts,2,function(x) sum(x[1:3])*2)
  sums_ref <- n - sums_alt
  freq_alt <- sums_alt/n
  freq_ref <- 1- freq_alt
#  het_count <- apply(counts,2,function(x) x[2]) # This should be simply a row from the counts matrix
  het_obs <- apply(counts,2,function(x) x[2]/sum(x[1:3]))
  #het_exp <- 2 * sums_alt/n * sums_ref/n
  return (list(#sums_alt = sums_alt,
               #sums_ref = sums_ref,
               freq_alt = freq_alt,
               freq_ref = freq_ref,
               #het_count = het_count,
               n = n,
               het_obs = het_obs))
}


pairwise_pop_fst_nei83 <- function(.x, by_locus = FALSE){
  # get the populations
  .group_levels = .x %>% group_keys()
  # create all combinations
  pairwise_combn <- utils::combn(nrow(.group_levels),2)
  # vector and matrix to store Fst for total and by locus
  Fst_tot <- rep(NA_real_, ncol(pairwise_combn))
  if (by_locus){
    Fst_locus <- matrix(NA_real_, nrow = count_loci(.x), ncol = ncol(pairwise_combn))
  }
  # summarise population frequencies
  pop_freqs_df <- group_map(.x, .f=~.gt_pop_freqs(.x))
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]

    n <-cbind(pop_freqs_df[[pop1]]$n,pop_freqs_df[[pop2]]$n)/2
    sHo <-cbind(pop_freqs_df[[pop1]]$het_obs,pop_freqs_df[[pop2]]$het_obs)
    mHo <- apply(sHo, 1, mean, na.rm = TRUE)
    freq_alt <- cbind(pop_freqs_df[[pop1]]$freq_alt, pop_freqs_df[[pop2]]$freq_alt)
    freq_ref <- cbind(pop_freqs_df[[pop1]]$freq_ref, pop_freqs_df[[pop2]]$freq_ref)

    # sum of squared frequencies
    sp2 <- freq_alt^2+freq_ref^2
    Hs <- (1 - sp2 - sHo/2/n)
    Hs <- n/(n - 1) * Hs
    np <- apply(n, 1, fun <- function(x) sum(!is.na(x)))
    # mean sample size over the populations
    mn <- apply(n, 1, fun <- function(x) {
      sum(!is.na(x))/sum(1/x[!is.na(x)])
    })
    # mean sum of square frequencies
    msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
    mp2 <- rowMeans(freq_alt)^2+rowMeans(freq_ref)^2
    mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
    Ht <- 1 - mp2 + mHs/mn/np - mHo/2/mn/np

    Dst <- Ht - mHs
    if (by_locus){
      Fst_locus[,i_col] = Dst/Ht
    }
    Fst_tot[i_col]<-mean(Dst)/mean(Ht)
  }
  # format nicely the objects
  group_combinations <- cbind(.group_levels[pairwise_combn[1,],],.group_levels[pairwise_combn[2,],])
  names(group_combinations) <- c(paste0(dplyr::group_vars(.x),"_1"),paste0(dplyr::group_vars(.x),"_2"))
  Fst_tot <- data.frame(group_combinations,value=Fst_tot)
  if (by_locus){
    rownames(Fst_locus)<-loci_names(.x)
    colnames(Fst_locus)<- apply(group_combinations,1,function(x)paste(x,collapse = "."))
  }
  if (by_locus){
    return(list(Fst_by_locus = Fst_locus, Fst_total = Fst_tot))
  } else{
    return(Fst_tot)
  }
}

pairwise_pop_fst_nei83(test_gt, by_locus = TRUE)



