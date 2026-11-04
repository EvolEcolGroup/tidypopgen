# to be sourced from within the vignette/articles directory
setwd(getSrcDirectory(function(x) {
  x
}))
unlink("./figure/sfs", recursive = TRUE)
knitr::knit("sfs.Rmd.orig", "sfs.Rmd")
