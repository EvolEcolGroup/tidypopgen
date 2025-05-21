# to be sourced from within the vignette/articles directory
setwd(getSrcDirectory(function(x) {
  x
}))
unlink("./figure/pca_aDNA_projection")
knitr::knit("pca_aDNA_projection.Rmd.orig", "pca_aDNA_projection.Rmd")
