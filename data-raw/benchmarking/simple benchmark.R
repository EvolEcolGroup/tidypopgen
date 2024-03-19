library(tictoc)
bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"))


tic()
missing_gt %>% loci_missingness()
toc()
