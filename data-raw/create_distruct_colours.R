## code to prepare `distruct_colours` dataset
distruct_colours <- read.csv(
  "./data-raw/files/default_colors.txt",
  header = FALSE
)$V1
# the distruct palette repeats 4 of its 60 colours (positions 35, 46, 54 and
# 60 repeat 4, 24, 12 and 59): keep each colour once, in the original order
distruct_colours <- unique(distruct_colours)
usethis::use_data(distruct_colours, overwrite = TRUE)
