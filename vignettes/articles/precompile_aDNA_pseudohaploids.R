# to be sourced from within the vignette/articles directory
setwd(getSrcDirectory(function(x) {
  x
}))
unlink("./figure/aDNA_pseudohaploids")
knitr::knit("aDNA_pseudohaploids.Rmd.orig", "aDNA_pseudohaploids.Rmd")
