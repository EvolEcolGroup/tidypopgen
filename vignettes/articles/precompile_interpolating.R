# to be sourced from within the vignette/articles directory
setwd(getSrcDirectory(function(x) {
  x
}))
unlink("./figure/interpolating", recursive = TRUE)
knitr::knit("interpolating.Rmd.orig", "interpolating.Rmd")
