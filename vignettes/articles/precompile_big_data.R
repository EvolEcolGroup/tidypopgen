# to be sourced from within the vignette/articles directory
setwd(getSrcDirectory(function(x) {
  x
}))
unlink("./figure/big_data", recursive = TRUE)
knitr::knit("big_data.Rmd.orig", "big_data.Rmd")
