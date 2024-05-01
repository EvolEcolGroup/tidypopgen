# internal functions that make life easier
.gt_bigsnp_cols <- function(.x){
  show_loci(.x)$big_index
}

.gt_bigsnp_rows <- function(.x){
  vctrs::vec_data(.x$genotypes)
}

.gt_get_bigsnp<-function(.x){
  attr(.x$genotypes,"bigsnp")
}


# a developer function to create various count summaries of a population, used to
# compute more complex statistics (e.g. pairwise fst, etc.).
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
