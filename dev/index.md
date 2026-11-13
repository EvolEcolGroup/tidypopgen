# tidypopgen

The goal of `tidypopgen` is to provide a tidy grammar of population
genetics, facilitating the manipulation and analysis of biallelic single
nucleotide polymorphisms (SNPs). `tidypopgen` scales to very large
genetic datasets by storing genotypes on disk, and performing operations
on them in chunks, without ever loading all data in memory. The full
functionalities of the package are described in [Carter et
al. (2025)](https://doi.org/10.1111/2041-210x.70204). Please cite this
paper if you use `tidypopgen` in your research.

## Installation

You can install the release version of `tidypopgen` from CRAN:

``` r

install.packages("tidypopgen")
```

You can install the latest development version directly from r-universe
(recommended):

``` r

install.packages('tidypopgen', repos = c('https://evolecolgroup.r-universe.dev',
                 'https://cloud.r-project.org'))
```

Alternatively, you can install `tidypopgen`using `devtools` (but you
might need to set up your development environment, which can be a bit
more complex):

``` r

install.packages("pak")
pak::pak("EvolEcolGroup/tidypopgen@dev")
```

## Examples

There are several vignettes designed to teach you how to use
`tidypopgen`. A short introduction to the package is available in the
[‘introduction’
vignette](https://evolecolgroup.github.io/tidypopgen/articles/tidypopgen.html).
A more detailed and technical description of the grammar of population
genetics, explaining how to manipulate individuals and loci, is
available in the [‘grammar’
vignette](https://evolecolgroup.github.io/tidypopgen/articles/a01_overview.html).

The [‘quality control’
vignette](https://evolecolgroup.github.io/tidypopgen/articles/a02_qc.html)
illustrates the `tidypopgen` functions that help running a full QC of a
dataset before analysis.

The [‘population genetic analysis’
vignette](https://evolecolgroup.github.io/tidypopgen/articles/a03_example_clustering_and_dapc.html)
provides a fully annotated example of how to run various population
genetics analyses with `tidypopgen`.

We also provide a [‘PLINK
cheatsheet’](https://evolecolgroup.github.io/tidypopgen/articles/a99_plink_cheatsheet.html)
aimed at translating common tasks performed in PLINK into `tidypopgen`
commands.

There is also an article showing how to manage aDNA samples that have
been coded as pseudohaploids, including how to project ancient DNA data
onto a PCA fitted to modern data and prepare data for admixtools: [‘aDNA
pseudohaploids’
article](https://evolecolgroup.github.io/tidypopgen/articles/aDNA_pseudohaploids.html).

Other articles cover specific analyses:

- the [‘spatial interpolation’
  article](https://evolecolgroup.github.io/tidypopgen/articles/interpolating.html)
  interpolates principal component scores across geographic space, using
  the integration of `tidypopgen` with `sf`;
- the [‘runs of homozygosity’
  article](https://evolecolgroup.github.io/tidypopgen/articles/roh_hgdp.html)
  finds and summarises the runs of homozygosity of the HGDP, with the
  same algorithm and results as PLINK 1.9 (`--homozyg`);
- the [‘site frequency spectra’
  article](https://evolecolgroup.github.io/tidypopgen/articles/sfs.html)
  computes single and joint site frequency spectra, dealing with missing
  genotypes by projection (as `easySFS` does), as the input of
  demographic inference with `dadi`, `moments` or `fastsimcoal2`.

Finally, `tidypopgen` is fast and can handle large datasets easily. See
a [‘benchmark’
article](https://evolecolgroup.github.io/tidypopgen/articles/benchmark_hgdp.html)
using the HGDP, a dataset of over 1000 individuals typed for 650k SNPs.
We can load the data, clean it, run imputation, PCA and pairwise Fst
among 51 populations in less than 20 seconds on a powerful desktop (and
less than a minute on a laptop). For even larger datasets, the [‘large
VCFs’
article](https://evolecolgroup.github.io/tidypopgen/articles/big_data.html)
merges many whole-genome VCF files into a single `gen_tibble`, storing
the names of the loci as integers to save memory, and the [‘large
RAD-seq datasets’
article](https://evolecolgroup.github.io/tidypopgen/articles/rad_seq.html)
handles de novo RAD-seq assemblies, in which every RAD locus is its own
chromosome, by storing the chromosomes as integers too.

## When something does not work

If something does not work, check the [issues on
GitHub](https://github.com/EvolEcolGroup/tidypopgen/issues) to see
whether the problem has already been reported. If not, feel free to
create an new issue. Please make sure you have updated to the latest
version of `tidypopgen` on r-universe/Github, as well as updating all
other packages on your system, and provide [a reproducible
example](https://reprex.tidyverse.org/) for the developers to
investigate the problem. Ideally, try to create a minimalistic dataset
that reproduces the error, as it will be much easier (and thus faster!)
for the developers to track down the problem.
