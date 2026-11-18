# A small VCF of mixed ploidy (Arabidopsis arenosa), the first 200 loci of the
# VCF of the polyploid tutorial of Bohutinska et al. (2023), see
# misc/AdegenetPolyplChapter_url.txt. The data are from Wos G. et al. (2019)
# Role of ploidy in colonization of alpine habitats in natural populations of
# Arabidopsis arenosa. Annals of Botany 124: 255-268. doi:10.1093/aob/mcz070
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6758580/).
ploidy_vcf <- vcfR::read.vcfR(
  "arenosa_snp_raw.fourfold.filtered.PolyplChapter.vcf.gz",
  nrow = 200
)
vcfR::write.vcf(ploidy_vcf, file = "arenosa_ploidy_example.vcf.gz")
