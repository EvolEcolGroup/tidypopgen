# tidypopgen <img src="./man/figures/logo.png" align="right" width="150"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/EvolEcolGroup/tidypopgen/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EvolEcolGroup/tidypopgen/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/EvolEcolGroup/tidypopgen/branch/main/graph/badge.svg?token=KLOzxJoLBO)](https://app.codecov.io/gh/EvolEcolGroup/tidypopgen)
<!-- badges: end -->

The goal of `tidypopgen` is to provide a tidy grammar of population genetics, facilitating 
the manipulation and analysis of biallelic single nucleotide 
polymorphisms (SNPs). `tidypopgen` scales to very large genetic datasets by storing 
genotypes on disk, and performing operations on them in chunks, without
ever loading all data in memory.

## Installation

You can install the latest version of `tidypopgen` directly from r-universe (recommended):
``` r
install.packages('tidypopgen', repos = "https://evolecolgroup.r-universe.dev")
```

Alternatively, you can install `tidypopgen`using `devtools` (but you might need to set up your development environment, 
which can be a bit more complex):
``` r
install.packages("devtools")
devtools::install_github("EvolEcolGroup/tidypopgen")
```

## Examples

There are a several vignettes designed to teach you how to use `tidypopgen`. 
The ['overview' vignette](https://evolecolgroup.github.io/tidypopgen/dev/articles/tidypopgen.html) 
explains how the data structures are designed and provides an illustration
of the grammar used to manipulate individuals and loci. 

The ['quality control' vignette](https://evolecolgroup.github.io/tidypopgen/dev/articles/a02_qc.html)
illustrates the `tidypopgen` functions that help
running a full QC of a dataset before analysis.

The ['population genetic analysis' vignette](https://evolecolgroup.github.io/tidypopgen/dev/articles/a03_example_clustering_and_dapc.html)
provides a fully annotated example of how to 
run various population genetics analyses with `tidypopgen`.

We also provide a ['PLINK cheatsheet'](https://evolecolgroup.github.io/tidypopgen/dev/articles/a99_plink_cheatsheet.html)
aimed at translating common tasks
performed in PLINK into `tidypopgen` commands.

Finally, `tidypopgen` is fast and can handle large dataset easily. See a 
[benchmark](https://evolecolgroup.github.io/tidypopgen/dev/articles/benchmark_hgdp.html) using the HGDP,
a dataset of over 1000 individuals typed for 650k SNPs. We can load the data, clean it,
run imputation, PCA and pairwise Fst among 51 populations in less than 20 seconds on a
powerful desktop.

## When something does not work

If something does not work, check the [issues on
GitHub](https://github.com/EvolEcolGroup/tidypopgen/issues) to see whether
the problem has already been reported. If not, feel free to create an
new issue. Please make sure you have updated to the latest version of
`tidypopgen` on r-universe/Github, as well as updating all other packages on your
system, and provide [a reproducible
example](https://reprex.tidyverse.org/)
for the developers to investigate the problem. Ideally, try to create a minimalistic
dataset that reproduces the error, as it will be much easier (and thus faster!)
for the developers to track down the problem. 
