# to be sourced from within the vignette/articles directory
setwd(getSrcDirectory(function(x) {
  x
}))
unlink("./figure/interpolating_hgdp")
knitr::knit("interpolating_hgdp.Rmd.orig", "interpolating_hgdp.Rmd")
