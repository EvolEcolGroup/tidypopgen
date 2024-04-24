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

