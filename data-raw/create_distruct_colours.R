## code to prepare `distruct_colours` dataset
distruct_colours <- read.csv("./data-raw/files/default_colors.txt",
                            header = FALSE)$V1
usethis::use_data(distruct_colours, overwrite = TRUE)
