# Compute Tajima's D with block jackknife standard errors

This function computes genome-wide Tajima's D, together with its
standard error estimated with a weighted block jackknife.

## Usage

``` r
jackknife_pop_tajimas_d(
  .x,
  type = c("matrix", "tidy", "list"),
  block_size,
  size_unit = c("snp", "cM", "bp")
)
```

## Arguments

- .x:

  a `gen_tibble` object, optionally grouped by one variable

- type:

  type of object to return, if using grouped method. One of "matrix",
  "tidy", or "list". Default is "matrix". Data can be grouped by more
  than one variable: groups are then labelled by the values of the
  grouping variables pasted with "." (e.g. `pop1.F`), and the "tidy"
  output has one column per grouping variable rather than a single
  `group` column. The `group` column, or the columns of the grouping
  variables, keep the type of the grouping variables (e.g. a factor
  keeps its levels), so that groups stay in their order when the output
  is sorted or plotted.

- block_size:

  The size of the blocks to leave out.

- size_unit:

  Either "snp", "cM" or "bp". If "snp", the block size is in number of
  SNPs. If "cM", the block size is in centiMorgans. If "bp", the block
  size is in base pairs.

## Value

if data is not grouped, a data frame with the following columns:

- `estimate`: Tajima's D computed on all loci

- `jack_estimate`: the bias corrected jackknife estimate

- `se`: the jackknife standard error

- `z`: the z score (jack_estimate/se)

- `n_blocks`: the number of blocks used if data are grouped, either: a
  data frame with the following columns:

- `quantity`: the quantity reported in the row (`estimate`,
  `jack_estimate`, `se`, `z` or `n_blocks`, as described above)

- `group`: the value of the quantity for the group (there will be as
  many of these columns as groups in the gen_tibble, and they will be
  named by the grouping levels) a tidy tibble with the following
  columns:

- `group`: the name of the group

- `estimate`, `jack_estimate`, `se`, `z`, `n_blocks`: as described above
  or a list of data frames, one per group, with the same columns as for
  ungrouped data.

## Details

The genome is divided into non-overlapping blocks, which never span more
than one chromosome, and Tajima's D is recomputed leaving out one block
at a time. The standard error is computed with the weighted delete-m
jackknife of Busing et al. (1999), where each block is weighted by the
number of loci it contains. Tajima's D is a function of sums across loci
(of the nucleotide diversity, the number of segregating sites, and
quantities that depend on the sample size of each locus), so leaving out
a block is achieved by subtracting its sums, without having to rescan
the genotypes. See
[`jackknife_stats_generic()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/jackknife_stats_generic.md)
for details.

Missing genotypes are accounted for by using the number of non-missing
allele copies at each locus as its sample size, following Bailey et al.
(2025) as implemented in `pixy`. See
[`pop_tajimas_d()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_tajimas_d.md)
for details. Loci with no genotyped individuals are ignored. The
variance of Tajima (1989) is computed for the mean sample size across
all loci (rounded to the nearest integer), and the same sample size is
used when leaving out each block: recomputing the rounded mean for each
block left out would make Tajima's D jump whenever this changes the
rounding (e.g. when the mean sample size is close to x.5), which the
jackknife would amplify into an inflated standard error and a biased
jackknife estimate.

Note that the variance of Tajima (1989) assumes that there is no
recombination among loci. Over a whole genome, Tajima's D is thus not
standardised: it is approximately proportional to \\\hat{\theta}\_\pi /
\hat{\theta}\_W - 1\\, with a factor that depends on the sample size
(about 5 for 10 allele copies, 4.2 for 20, and 3.4 for 100), and with
genome-wide data its standard error is much smaller than 1. Use the
jackknife standard error, rather than the rules of thumb used for single
loci, to judge whether genome-wide Tajima's D differs from 0, and
compare groups with different sample sizes with caution. The z score
tests whether genome-wide Tajima's D differs from 0, which most
demographic histories other than a constant population size will cause.
Blocks should be larger than the extent of linkage disequilibrium, or
the standard error will be underestimated.

## References

Busing FMTA, Meijer E, Van Der Leeden R (1999) Delete-m jackknife for
unequal m. Statistics and Computing 9:3-8.

Tajima F (1989) Statistical method for testing the neutral mutation
hypothesis by DNA polymorphism. Genetics 123:585-595.

Bailey N, Stevison L, Samuk K (2025) Correcting for bias in estimates of
\\\theta_W\\ and Tajima's D from missing data in next-generation
sequencing. Molecular Ecology Resources, e14104.
[doi:10.1111/1755-0998.14104](https://doi.org/10.1111/1755-0998.14104)

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Calculate Tajima's D using blocks of 2 SNPs
example_gt %>%
  jackknife_pop_tajimas_d(
    block_size = 2,
    size_unit = "snp"
  )
#> # A tibble: 5 × 4
#>   quantity       pop1     pop2  pop3
#>   <chr>         <dbl>    <dbl> <dbl>
#> 1 estimate      0.944  -0.780   -Inf
#> 2 jack_estimate 1.04   -0.861     NA
#> 3 se            0.329   0.0297    NA
#> 4 z             3.14  -29.0       NA
#> 5 n_blocks      3       3          3
```
