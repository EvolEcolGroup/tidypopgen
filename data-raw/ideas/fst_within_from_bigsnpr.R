
# this is equivalent to --fst --within
# TODO this needs to be formatted in a manner similar to Hudson fst
pairwise_pop_fst_wc <- function(.x, by_locus=FALSE){
  # simple function to create a dataframe with af and N (allele freq and number of valid alleles)
  make_ac_df <- function(x){
    n_tot <- nrow(x)*2
    x %>% reframe(af = loci_maf(x),N= n_tot - loci_missingness(x, as_counts=TRUE)*2)
  }
  ac_list <- group_map(.x, .f=~make_ac_df(.x))
  bigsnpr::snp_fst(ac_list,overall=!by_locus)

}


