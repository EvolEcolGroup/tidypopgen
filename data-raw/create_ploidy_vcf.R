#File from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6758580/
ploidy_vcf <- vcfR::read.vcfR(
  "arenosa_snp_raw.fourfold.filtered.PolyplChapter.vcf.gz",
  nrow = 200
)
vcfR::write.vcf(ploidy_vcf, file = "arenosa_ploidy_example.vcf.gz")
