# to be sourced from within the vignette/articles directory; set the
# environment variable PLINK to the path of PLINK 1.9 (or have `plink` on the
# PATH) to compare with PLINK
setwd(getSrcDirectory(function(x) {
  x
}))
unlink("./figure/roh_hgdp", recursive = TRUE)
knitr::knit("roh_hgdp.Rmd.orig", "roh_hgdp.Rmd")
