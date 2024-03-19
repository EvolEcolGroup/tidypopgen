library(tictoc)
bed_file <- "./data-raw/benchmarking/hgdp_100k.bed"
hgdp_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"))
gt_save(hgdp_gt)

hdgp_sub <- hgdp_gt %>% select_loci(1:5000)
tic()
foo <- hgdp_gt %>% loci_hwe
toc()

tic()
foo2 <- hdgp_sub %>% loci_hwe()
toc()
