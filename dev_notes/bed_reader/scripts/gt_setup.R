# Make the gen_tibble of the baboon subset (new .gt format), and a copy of its
# .gt in the old format (gzfile level 1, character alleles)
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(tidypopgen, lib.loc = args[1]))
x <- gen_tibble(args[2], backingfile = args[3], names_as_int = TRUE,
                n_cores = 8, quiet = TRUE)
con <- gzfile(paste0(args[3], "_old.gt"), open = "wb", compression = 1)
saveRDS(x, con)
close(con)
# the two files give identical gen_tibbles
a <- gt_load(paste0(args[3], ".gt"))
b <- gt_load(paste0(args[3], "_old.gt"))
cat("identical loci:", identical(show_loci(a), show_loci(b)),
    "identical ids:", identical(a$id, b$id), "\n")
cat("sizes (GB): new", file.size(paste0(args[3], ".gt")) / 1e9,
    "old", file.size(paste0(args[3], "_old.gt")) / 1e9, "\n")
