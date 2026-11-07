# Detect runs of homozygosity with `detectRUNS` (deprecated)

This function is deprecated: use
[`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_indiv_roh.md),
which finds runs of homozygosity as PLINK (`--homozyg`), with its
parameters. It is kept so that old code still gives the same results: it
calls the sliding-window algorithm of the package `detectRUNS`
([`detectRUNS::slidingRUNS.run()`](https://rdrr.io/pkg/detectRUNS/man/slidingRUNS.run.html))
on each individual, as
[`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_indiv_roh.md)
used to do. It needs `detectRUNS` to be installed, and will be removed
in a future version of `tidypopgen`.

## Usage

``` r
gt_roh_window(
  .x,
  window_size = 15,
  threshold = 0.05,
  min_snp = 3,
  heterozygosity = FALSE,
  max_opp_window = 1,
  max_miss_window = 1,
  max_gap = 10^6,
  min_length_bps = 1000,
  min_density = 1/1000,
  max_opp_run = NULL,
  max_miss_run = NULL
)
```

## Arguments

- .x:

  a
  [gen_tibble](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)

- window_size:

  the size of sliding window (number of SNP loci) (default = 15)

- threshold:

  the threshold of overlapping windows of the same state
  (homozygous/heterozygous) to call a SNP in a RUN (default = 0.05)

- min_snp:

  minimum n. of SNP in a RUN (default = 3)

- heterozygosity:

  should we look for runs of heterozygosity (instead of homozygosity?
  (default = FALSE)

- max_opp_window:

  max n. of SNPs of the opposite type (e.g. heterozygous snps for runs
  of homozygosity) in the sliding window (default = 1)

- max_miss_window:

  max. n. of missing SNP in the sliding window (default = 1)

- max_gap:

  max distance between consecutive SNP to be still considered a
  potential run (default = 10^6 bps)

- min_length_bps:

  minimum length of run in bps (defaults to 1000 bps = 1 kbps)

- min_density:

  minimum n. of SNP per kbps (defaults to 0.1 = 1 SNP every 10 kbps)

- max_opp_run:

  max n. of opposite genotype SNPs in the run (optional)

- max_miss_run:

  max n. of missing SNPs in the run (optional)

## Value

A dataframe with RUNs of Homozygosity or Heterozygosity in the analysed
dataset, in the format of `detectRUNS`, with columns "group", "id",
"chrom", "nSNP", "from", "to", "lengthBps" (group: population, breed,
case/control etc.; id: individual identifier; chrom: chromosome on which
the run is located; nSNP: number of SNPs in the run; from: starting
position of the run, in bps; to: end position of the run, in bps;
lengthBps: size of the run)

## Details

If the
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
is grouped, then the grouping variable is used to fill in the 'group'
column (with several grouping variables, their values are pasted with
".", e.g. `pop1.F`). Otherwise, the 'group' column is filled with the
same values as the 'id' column.

## References

Biscarini F, Cozzi P, Gaspa G, Marras G (2019) `detectRUNS`: Detect Runs
of Homozygosity and Runs of Heterozygosity in Diploid Genomes. R
package.
[doi:10.32614/CRAN.package.detectRUNS](https://doi.org/10.32614/CRAN.package.detectRUNS)

## See also

[`detectRUNS::slidingRUNS.run()`](https://rdrr.io/pkg/detectRUNS/man/slidingRUNS.run.html)
which this function wraps.

## Examples

``` r
bed_file <- system.file("extdata", "roh", "roh_sim.bed",
  package = "tidypopgen"
)
sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
# deprecated: use windows_indiv_roh() instead
sim_roh <- suppressWarnings(gt_roh_window(sim_gt %>% group_by(population)))
#> No RUNs found for animal outbred_01
#> No RUNs found for animal outbred_02
#> N. of RUNS for individual outbred_03 is: 1
#> No RUNs found for animal outbred_04
#> N. of RUNS for individual outbred_05 is: 1
#> No RUNs found for animal outbred_06
#> N. of RUNS for individual outbred_07 is: 3
#> No RUNs found for animal outbred_08
#> N. of RUNS for individual outbred_09 is: 2
#> N. of RUNS for individual outbred_10 is: 2
#> N. of RUNS for individual moderate_01 is: 3
#> N. of RUNS for individual moderate_02 is: 3
#> N. of RUNS for individual moderate_03 is: 3
#> N. of RUNS for individual moderate_04 is: 1
#> N. of RUNS for individual moderate_05 is: 2
#> N. of RUNS for individual moderate_06 is: 1
#> N. of RUNS for individual moderate_07 is: 2
#> N. of RUNS for individual moderate_08 is: 3
#> N. of RUNS for individual moderate_09 is: 3
#> N. of RUNS for individual moderate_10 is: 2
#> N. of RUNS for individual inbred_01 is: 4
#> N. of RUNS for individual inbred_02 is: 3
#> N. of RUNS for individual inbred_03 is: 1
#> N. of RUNS for individual inbred_04 is: 4
#> N. of RUNS for individual inbred_05 is: 5
#> N. of RUNS for individual inbred_06 is: 4
#> N. of RUNS for individual inbred_07 is: 3
#> N. of RUNS for individual inbred_08 is: 3
#> N. of RUNS for individual inbred_09 is: 7
#> N. of RUNS for individual inbred_10 is: 4
head(sim_roh)
#>     group         id chrom nSNP    from      to lengthBps
#> 1 outbred outbred_03     3   16 1546240 1628083     81843
#> 2 outbred outbred_05     1   15  529488  686908    157420
#> 3 outbred outbred_07     1   20 1637458 1810256    172798
#> 4 outbred outbred_07     2   15  956322 1029451     73129
#> 5 outbred outbred_07     2   15 5651823 5751851    100028
#> 6 outbred outbred_09     1   16 7929001 8109667    180666
```
