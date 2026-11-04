# Create the lists of SGDP individuals used in the big data article
# (vignettes/articles/big_data.Rmd.orig): the first 10 fully public
# individuals from each of five regions. Run this script from the root of the
# repository, BEFORE `s02_create_sgdp_subset.sh`, which downloads the
# individuals listed in the `<region>_subset.txt` files written here, next to
# both scripts.

sgdp_dir <- "data-raw/sgdp"
metadata <-
  read.table(
    file.path(sgdp_dir, "SGDP_metadata.279public.21signedLetter.44Fan.samples.csv"),
    header = TRUE,
    sep = ";"
  )

# subset to public data
metadata <- subset(metadata, Embargo == "FullyPublic")

# the regions to sample, named by the prefix of their file
regions <- c(
  afr = "Africa",
  west_eur = "WestEurasia",
  america = "America",
  oceania = "Oceania",
  east_asia = "EastAsia"
)
for (prefix in names(regions)) {
  region_ids <- metadata$Illumina_ID[metadata$Region == regions[[prefix]]]
  write.table(region_ids[1:10],
    file.path(sgdp_dir, paste0(prefix, "_subset.txt")),
    col.names = FALSE, quote = FALSE, row.names = FALSE
  )
}
