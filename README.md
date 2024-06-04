# tidypopgen

<!-- badges: start -->
  [![R-CMD-check](https://github.com/EvolEcolGroup/tidypopgen/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EvolEcolGroup/tidypopgen/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

The goal of `tidypopgen` is to provide a tidy grammar of population genetics, facilitating 
the manipulation and analysis of genetic data. Currently, it is focussed on biallelic single nucleotide
polymorphisms (SNPs).

## Installation

You need `devtools` to install `tidypopgen`. If you haven't done so already, install it with:
``` r
install.packages("devtools")
```

You are then ready to install the development version of `tidypopgen` from [GitHub](https://github.com/) with:
``` r
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

