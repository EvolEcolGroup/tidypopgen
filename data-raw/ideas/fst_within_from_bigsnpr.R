
# this is equivalent to --fst --within
# TODO this needs to be formatted in a manner similar to Hudson fst
pairwise_pop_fst_wc <- function(.x, by_locus=FALSE){
  # simple function to create a dataframe with af and N (allele freq and number of valid alleles)
  make_ac_df <- function(x){
    n_tot <- nrow(x)*2
    x %>% reframe(af = loci_maf(x),N= n_tot - loci_missingness(x, as_counts=TRUE)*2)
  }
  ac_list <- group_map(.x, .f=~make_ac_df(.x))
  browser()
  bigsnpr::snp_fst(ac_list,overall=!by_locus)

}


x <- test_gt %>% filter(population=="pop1")

make_ac_df <- function(.x){
  counts <- bigstatsr::big_counts( .gt_get_bigsnp(.x)$genotypes,
                                   ind.row =.gt_bigsnp_rows(.x),
                                   ind.col = .gt_bigsnp_rows(.x))
  sums <- apply(counts,2,function(x) x[2]+2*x[3])
  n <- apply(counts,2,function(x) sum(x[1:3])*2)
  het_obs <- apply(counts,2,function(x) x[2]/sum(x[1:3]))
  het_exp <- 2 * sums/n * (1-sums/n)
}
