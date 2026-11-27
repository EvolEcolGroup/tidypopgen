# to be sourced from within the vignette/articles directory
setwd(getSrcDirectory(function(x) {
  x
}))
unlink("./figure/rad_seq", recursive = TRUE)
knitr::knit("rad_seq.Rmd.orig", "rad_seq.Rmd")
