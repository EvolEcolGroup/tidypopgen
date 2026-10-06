metadata <-
  read.table("data-raw/SGDP_metadata.279public.21signedLetter.44Fan.samples.csv",
             header = TRUE,
             sep = ";")

# subset to public data
metadata <- subset(metadata, Embargo == "FullyPublic")
# subset to African
metadata_afr <- subset(metadata, Region == "Africa")
afr_subsample <- metadata_afr["Illumina_ID"][c(1:10),]
write.table(afr_subsample, col.names = FALSE, quote = FALSE, row.names = FALSE, "data-raw/afr_subset.txt")


# subset to European/West Eurasian
metadata_west_eur <- subset(metadata, Region == "WestEurasia")
west_eur_subsample <- metadata_west_eur["Illumina_ID"][c(1:10),]
write.table(west_eur_subsample, col.names = FALSE, quote = FALSE, row.names = FALSE, "data-raw/west_eur_subset.txt")

# subset to America
metadata_america <- subset(metadata, Region == "America")
america_subsample <- metadata_america["Illumina_ID"][c(1:10),]
write.table(america_subsample,  col.names = FALSE, quote = FALSE, row.names = FALSE, "data-raw/america_subset.txt")


# subset to Oceania
metadata_oceania <- subset(metadata, Region == "Oceania")
oceania_subsample <- metadata_oceania["Illumina_ID"][c(1:10),]
write.table(oceania_subsample, col.names = FALSE,  quote = FALSE, row.names = FALSE, "data-raw/oceania_subset.txt")



