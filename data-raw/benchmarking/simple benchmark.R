library(tictoc)
bed_file <- "./data-raw/benchmarking/hgdp_100k.bed"
hgdp_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"))
hgdp_gt <- gt_impute_simple(hgdp_gt)
tic()
test_ibs <- snp_ibs(gt_get_bigsnp(hgdp_gt)$genotypes, block.size = 10000)
toc()

hdgp_sub <- hgdp_gt %>% select_loci(1:5000)
tic()
foo <- hgdp_gt %>% loci_hwe
toc()

tic()
foo2 <- hdgp_sub %>% loci_hwe()
toc()


foo <- snp_readBed2("./data-raw/benchmarking/hgdp_100k.bed", backingfile = tempfile("missing_"),
                    ind.col = 1:30000,
                    ncores = 1
)
test_snp <-big_attach(foo)
library(tictoc)
tic()
test_ibs <- snp_ibs(test$genotypes)
toc()
