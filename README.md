# tidypopgen <img src="./man/figures/logo.png" align="right" width="150"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/EvolEcolGroup/tidypopgen/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EvolEcolGroup/tidypopgen/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/EvolEcolGroup/tidypopgen/branch/main/graph/badge.svg?token=KLOzxJoLBO)](https://app.codecov.io/gh/EvolEcolGroup/tidypopgen)
<!-- badges: end -->

The goal of `tidypopgen` is to provide a tidy grammar of population genetics, facilitating 
the manipulation and analysis of genetic data. Currently, it is focussed on biallelic single nucleotide
polymorphisms (SNPs).

We are making available a **preview** version of the package. Everything should work,
but be vigilant as not all functions have been tested extensively. Also, the syntax of
functions is still slightly in flux, so the names and order of some parameters might
change. Make sure you don't rely on position when giving multiple parameters to
a function, so that any change will trigger an informative error in the code that 
will allow you to update.

## Installation

You can install the latest version of `tidypopgen` directly from r-universe (reccomended):
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
The
'overview' vignette explains how the data structures are designed, and provides a illustration
of the grammar used to manipulate individuals and loci. 

The 'example workflow' vignette provides a fully annotated example of how to 
run population genetics analysis with `tidypopgen`.

The 'quality control' vignette illustrates the `tidypopgen` functions that help
running a fully QC of a dataset before analysis.

Finally, we provide a 'PLINK cheatsheet' aimed at translating common tasks
performed in PLINK into `tidypopgen` commands.

## When something does not work

This is a **preview** version of `tidypopgen`. Everything should work, but not all
elements have been extensively tested. If something does not work, check the [issues on
GitHub](https://github.com/EvolEcolGroup/tidypopgen/issues) to see whether
the problem has already been reported. If not, feel free to create an
new issue. Please make sure you have updated to the latest version of
`tidypopgen` on Github, as well as updating all other packages on your
system, and provide [a reproducible
example](https://reprex.tidyverse.org/)
for the developers to investigate the problem. Ideally, try to create a minimalistic
dataset that reproduces the error, as it will be much easier (and thus faster!)
for the developers to track down the problem. 
