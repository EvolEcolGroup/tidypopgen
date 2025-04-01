# to be sourced from within the vignette/articles directory
setwd(getSrcDirectory(function(x) {x}))
unlink("./figure")
knitr::knit("benchmark_hgdp.Rmd.orig", "benchmark_hgdp.Rmd")
