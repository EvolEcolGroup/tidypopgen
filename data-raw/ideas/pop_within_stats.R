.x <- pinus_gt
.x <- .x[1:5,] %>% select_loci(1:10)
# we compute quantities by locus (the can then be combined for individuals later)
# counts of genotypes
.counts <- bigstatsr::big_counts(.gt_get_bigsnp(.x)$genotypes,
                                 ind.row = .gt_bigsnp_rows(.x),
                                 ind.col = .gt_bigsnp_cols(.x))
# number of alleles
.n_k <- rbind(.counts[1,]*2+.counts[2,],
              .counts[3,]*2+.counts[2,])
# remove invariant sites
to_remove <- (.n_k[1,]==0 | .n_k[2,]==0)
.counts <- .counts[,!to_remove]
.n_k <- .n_k[,!to_remove]

# total sample size
.N <- apply(.counts[1:3,],2,sum)
# frequencies of alleles
.p <- .n_k[1,]/(.N*2)
.q <- 1- .p
# observed heterozygosity
h_o <- .counts[2,]/.N
# unbiased within pop gene diversity (h_s)
h_s <- (1-(.p^2+.q^2))*(2*.N/(2*.N-1))

Fs <- (h_s - h_o)/h_s
# if it is a grouped tibble with multiple populations

# overall gene diversity



##########################################
# convenient functs
.gt_bigsnp_cols <- function(.x){
  show_loci(.x)$big_index
}

.gt_bigsnp_rows <- function(.x){
  vctrs::vec_data(.x$genotypes)
}
