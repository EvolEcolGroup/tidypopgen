# set the working directory to the current script
setwd(getSrcDirectory(function(x) {
  x
}))
# move existing figures into a backup folder
fs::dir_copy("./figure",
             "./figure_backup",
             overwrite = TRUE)
rmd_orig <- readLines("./benchmark_hgdp.Rmd.orig")
# add text about the laptop
rmd_orig <- append(rmd_orig,
        values = c("# Running the benchmark on a laptop",
                   "We will now run this benchmark on a laptop with a `r benchmarkme::get_cpu()` CPU",
                   "and `r as.character(print(benchmarkme::get_ram(), unit_system = 'iec'))` of RAM,",
                   "limiting the number of cores to 4."),
       after = grep("# Summary", rmd_orig)+1)
# replace setting number of threads to 4
rmd_orig <- gsub("n_cores <- 20", "n_cores <- 4", rmd_orig)
# write out the file
writeLines(rmd_orig, "./laptop_benchmark_hgdp.Rmd.orig")
# run the benchmark
knitr::knit("laptop_benchmark_hgdp.Rmd.orig", "laptop_benchmark_hgdp.Rmd")
# ditch the new figures and bring back to old ones
unlink("./figure")
fs::dir_copy("./figure_backup",
             "./figure",
             overwrite = TRUE)
unlink("./figure_backup", recursive = TRUE)
# append the new benchmark values to the benchmark file
rmd_bench <- readLines("./benchmark_hgdp.Rmd")
rmd_laptop <- readLines("./laptop_benchmark_hgdp.Rmd")
# subset to the last summary section
rmd_laptop <- rmd_laptop[grep("# Running the benchmark on a laptop", rmd_laptop):length(rmd_laptop)]
rmd_bench <- c(rmd_bench, rmd_laptop)
# write out the file
writeLines(rmd_bench, "./benchmark_hgdp.Rmd")
# and now remove the laptop benchmark files
unlink("./laptop_benchmark_hgdp.Rmd")
unlink("./laptop_benchmark_hgdp.Rmd.orig")
