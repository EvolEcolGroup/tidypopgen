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
                        genetic_dist = as.integer(rep(0,6)),
                        allele_ref = c("A","T","C","G","C","T"),
                        allele_alt = c("T","C", NA,"C","G","A"))

test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)
test_gt <- test_gt %>% group_by(population)

# convert to hierfstat
test_hier <- gt_as_hierfstat(test_gt)

hier_basic <- hierfstat::basic.stats(test_hier)
# note that Fstp, FIS and Dest are not simply averages

hier_fst_wc <- hierfstat::pairwise.WCfst(test_hier)
hier_fst_nei <- hierfstat::pairwise.neifst(test_hier)

# a developer function to create various count summaries of a population, used to
# compute more complex statistics (e.g. pairwise fst, etc.).
# Specifically, we compute:
# sums_alt sum alleles for alt
# sums_ref sums alleles for ref
# n sums of all alleles
# het_obs observed heterozygosity
# het_exp expected heterozygosity
.gt_pop_counts <- function(.x){
  counts <- bigstatsr::big_counts( .gt_get_bigsnp(.x)$genotypes,
                                   ind.row =.gt_bigsnp_rows(.x),
                                   ind.col = .gt_bigsnp_cols(.x))
  sums_alt <- apply(counts,2,function(x) x[2]+2*x[3])
  n <- apply(counts,2,function(x) sum(x[1:3])*2)
  sums_ref <- n - sums_alt
  het_obs <- apply(counts,2,function(x) x[2]/sum(x[1:3]))
  het_exp <- 2 * sums_alt/n * sums_ref/n
  return (list(sums_alt = sums_alt,
               sums_ref = sums_ref,
               n = n,
               het_obs = het_obs,
               het_exp = het_exp))
}

pairwise_pop_fst_nei73 <- function(.x, by_locus = FALSE){
  browser()
  pop_counts_df <- group_map(.x, .f=~.gt_pop_counts(.x))
  # get the grouping column, and create all pairwise combination of indices
  .group_levels = .x %>% group_keys()
  pairwise_combn <- utils::combn(nrow(.group_levels),2)
  n_loci <- count_loci(.x)
  Hs_pair <- matrix(NA_real_, nrow = n_loci, ncol = ncol(pairwise_combn))
  Ht_pair <- matrix(NA_real_, nrow = n_loci, ncol = ncol(pairwise_combn))
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]
    Hs_pair[,i_col] <-  (pop_counts_df[[pop1]]$het_exp + pop_counts_df[[pop2]]$het_exp)/2
    # Ht_pair is 1-p_bar^2 - q_bar^2
    Ht_pair[, i_col] <-
      1 - ((pop_counts_df[[pop1]]$sums_alt + pop_counts_df[[pop2]]$sums_alt) /
             (pop_counts_df[[pop1]]$n + pop_counts_df[[pop2]]$n)) ^ 2 -
      ((pop_counts_df[[pop1]]$sums_ref + pop_counts_df[[pop2]]$sums_ref) /
         (pop_counts_df[[pop1]]$n + pop_counts_df[[pop2]]$n)) ^ 2
  }
  if (by_locus){
    1-(Hs_pair/Ht_pair)
    # tidy it properly
    #TODO
  } else {
    1-(colSums(Hs_pair,na.rm=TRUE)/colSums(Ht_pair,na.rm=TRUE))
    # TODO TIDY IT PROPERLY
  }
}

pairwise_pop_fst_wc(test_gt)
