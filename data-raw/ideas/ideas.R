# create a method is_diploid working on SNPbin and on list of SNPbin
# many statistics only make sense on diploid organisms
# hierfstat has very clear description of statistics on dosage data
# with explicit formulae

# Also look at package fst4pg
# And paper Estimating and interpresting Fst:The impact of rare variants advocating
# Hudson Fst

# Think about broom tidiers for various analysis
# https://www.tidymodels.org/learn/develop/broom/

# Look at
# https://www.tidymodels.org/learn/statistics/k-means/
# for a good example on k-means clustering (which is available from adegenet)

# For PCA, look at:
# https://broom.tidymodels.org/reference/tidy.prcomp.html
# https://clauswilke.com/blog/2020/09/07/pca-tidyverse-style/
# and https://cran.r-project.org/web/packages/ggbiplot/ggbiplot.pdf
# for plotting dudi.pca
# And some ideas from factorextra (not fully tidy, but some interesting ideas)
# http://www.sthda.com/english/wiki/get-pca-extract-the-results-for-individuals-variables-in-principal-component-analysis-r-software-and-data-mining

## A tutorial for analysis of SNP with adegenet that we could try to repeat with our pipeline:
#https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html

# Pipeline with adegenet and hierfstat to test pop gen summaries:
#https://popgen.nescent.org/DifferentiationSNP.html

# additional notes on stats
# https://www.uwyo.edu/dbmcd/popecol/maylects/fst.html

# reading vcf in chunks
#https://rdrr.io/bioc/genotypeeval/man/ReadVCFDataChunk.html

# https://github.com/Zilong-Li/vcfpp

##################################
#big SVD
# https://stackoverflow.com/questions/46253537/computing-the-null-space-of-a-bigmatrix-in-r/46380540#46380540

# really good tutorial on the maths of gwas
#https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS5.html

# Note  on running admixture and pca by John Novembre:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8722024/
