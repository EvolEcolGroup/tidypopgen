# to be sourced from within the vignette/articles directory
setwd(getSrcDirectory(function(x) {
  x
}))
unlink("./figure/pca_aDNA_projection")
unlink("./data/NearEastPublic", recursive = "TRUE")
unlink("./data/NearEastPublic.tar.gz")
knitr::knit("pca_aDNA_projection.Rmd.orig", "pca_aDNA_projection.Rmd")
