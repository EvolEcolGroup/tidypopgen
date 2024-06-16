# example script from
# https://github.com/mbohutinska/PolyplChapter/blob/main/PartB/AdegenetPolyplChapter.R

#set working directory - to your species folder (alternatively: RStudio -> Session -> Set Working Directory -> To Source File location)
setwd("/home/aa/popgenPolypl/")

#load libraries
library(adegenet)
library(vcfR)
library(gplots)
library(ggplot2)
library(StAMPP)
source("adegenet_functions.r")

############### SECTION 1 - GENETIC STRUCTURE ###############
#Here, you will infer the genetic structure among your populations, calculating their genetic distances. You will use functions from the package vcfR to handle the data and package adegenet and StAMPP to calculate the genetic distances and visualize


### A. import SNP data from VCF - necessary to change path using setwd()
vcf <- read.vcfR("arenosa_snp_raw.fourfold.filtered.PolyplChapter.vcf.gz")
aa.genlight <- vcfR2genlight.tetra(vcf)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_") # add real SNP.names
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3)  # add pop names



#Check if everything fits
indNames(aa.genlight)
pop(aa.genlight)
aa.genlight.noOut <- aa.genlight[pop=c(1,3:5)]
popNames(aa.genlight.noOut)
### B. calculate PCA

pca.1 <- glPcaFast(aa.genlight.noOut, nf=300)

#You can save it using the pdf function
pdf ("PCA_PolyplChapter_ax12.pdf", width=4, height=4)
plot(pca.1$scores, type="n", xlab=paste(round((pca.1$eig[1]/sum(pca.1$eig)*100),1)," %"), ylab=paste(round((pca.1$eig[2]/sum(pca.1$eig)*100),1)," %"))
s.class(pca.1$scores, pop(aa.genlight.noOut), add.plot = T, xax=1, yax=2,cstar = F,cellipse = F,cpoint = 4,col = c("blue3","red","blue","red3")) # insert your PCA plotting code with better color scheme here ###
dev.off()

### C. calculate genetic distances and visualize
aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)   # Nei's 1972 genetic distance between pops

#plot tree - neighbor joining tree
pdf ("njTree_PolyplChapter.pdf", width=4, height=4)
plot(nj(aa.D.pop))
dev.off()

# plot heatmap of the population distance matrix
colnames(aa.D.pop)<-rownames(aa.D.pop) # name the rows of a matrix

pdf("stamppdist.heatmap.PolypChapter.pdf", width=10,height=10)
heatmap.2(aa.D.pop, trace="none", cexRow=0.7, cexCol=0.7)
dev.off()

#Does the analysis of genetic structure suggest a single or parallel origin of tetraploids in your dataset?
