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
x <- test_gt %>% filter(population=="pop1")

make_ac_df <- function(.x){
  counts <- bigstatsr::big_counts( .gt_get_bigsnp(.x)$genotypes,
                                   ind.row =.gt_bigsnp_rows(.x),
                                   ind.col = .gt_bigsnp_cols(.x))
  sums <- apply(counts,2,function(x) x[2]+2*x[3])
  n <- apply(counts,2,function(x) sum(x[1:3])*2)
  het_obs <- apply(counts,2,function(x) x[2]/sum(x[1:3]))
  het_exp <- 2 * sums/n * (1-sums/n)
  return (list(sums = sums,
               n = n,
               het_obs = het_obs,
               het_exp = het_exp))
}

pairwise_pop_fst_wc <- function(.x){
  ac_df <- group_map(.x, .f=~make_ac_df(.x))
  # get the grouping column, and create all pairwise combination of indices
  .group_levels = .x %>% group_keys()
  pairwise_combn <- utils::combn(nrow(.group_levels),2)
  numerator <- matrix(NA_real_, nrow = nrow(pairwise_combn), ncol = n_loci)
  denominator <- matrix(NA_real_, nrow = nrow(pairwise_combn), ncol = n_loci)
  for (i_col in seq_len(ncol(pairwise_combn))){
    pop1 <- pairwise_combn[1,i_col]
    pop2 <- pairwise_combn[2,i_col]



  }
  browser()
}

pairwise_pop_fst_wc(test_gt)
