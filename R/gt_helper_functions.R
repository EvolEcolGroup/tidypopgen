# internal functions that make life easier
.gt_bigsnp_cols <- function(.x) {
  show_loci(.x)$big_index
}

.gt_bigsnp_rows <- function(.x) {
  vctrs::vec_data(.x$genotypes)
}

.gt_get_bigsnp <- function(.x) {
  # if this is a gen_tibble
  if (inherits(.x, "gen_tbl")) {
    attr(.x$genotypes, "bigsnp")
  } else {
    attr(.x, "bigsnp")
  }
}
