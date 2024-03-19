library(tictoc)
bed_file <- "./data-raw/benchmarking/hgdp_100k.bed"
hgdp_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"))
gt_save(hgdp_gt)

tic()
foo <- hgdp_gt %>% loci_missingness()
toc()

tic()
foo2 <- hgdp_gt %>% loci_missingness()
toc()
