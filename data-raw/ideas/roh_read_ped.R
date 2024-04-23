genotypeFilePath <- system.file(
  "extdata", "Kijas2016_Sheep_subset.ped", package="detectRUNS")
mapFilePath <- system.file(
  "extdata", "Kijas2016_Sheep_subset.map", package="detectRUNS")

sheep_ped <- snpStats::read.pedfile(genotypeFilePath)
map <- read.table(mapFilePath)
names(map)<-c("chromosome", "id", "genetic.distance","position")
sheep_ped$map<-cbind(sheep_ped$map,map)

snpStats::write.plink("./data-raw/datasets/sheep",snps = sheep_ped$genotypes,
                      subject.data = sheep_ped$fam,
                      chromosome = sheep_ped$map$chromosome,
                      allele.1 = sheep_ped$map$allele.1,
                      allele.2 = sheep_ped$map$allele.2,
                      position = sheep_ped$map$position,
                      genetic.distance = sheep_ped$map$genetic.distance,
                      human.genome = FALSE
                      )
bed_path <- "./data-raw/datasets/sheep.bed"
sheep_gt <- tidypopgen::gen_tibble(bed_path, backingfile = tempfile())
