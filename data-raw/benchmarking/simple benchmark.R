library(tictoc)
bed_file <- "./data-raw/benchmarking/hgdp_100k.bed"
hgdp_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"))
hgdp_gt <- gt_impute_simple(hgdp_gt)
tic()
test_ibs <- snp_ibs(gt_get_bigsnp(hgdp_gt)$genotypes, block.size = 10000)
toc()

tic()
test_ibs <- snp_ibs_loop(test,rows_along(test),cols_along(test),1)
toc()


hdgp_sub <- hgdp_gt %>% select_loci(1:5000)
tic()
foo <- hgdp_gt %>% loci_hwe
toc()

tic()
foo2 <- hdgp_sub %>% loci_hwe()
toc()


foo <- snp_readBed("./data-raw/benchmarking/hgdp_100k.bed", backingfile = tempfile("missing_"))
test_snp <-bigsnpr::snp_attachExtdata()
library(tictoc)
tic()
test_ibs <- snp_ibs(test_snp$genotypes)
toc()

tic()
test_ibs_R <- snp_ibs(test_snp$genotypes)
toc()


