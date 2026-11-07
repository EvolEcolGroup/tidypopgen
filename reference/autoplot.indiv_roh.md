# Autoplots for `indiv_roh` objects

Plots of runs of homozygosity found by
[`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_indiv_roh.md),
equivalent to the plots of `detectRUNS`. Runs are coloured (or split) by
group if the
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
was grouped. The following types of plots are available:

- `runs`: the runs of each individual along each chromosome (as
  [`detectRUNS::plot_Runs()`](https://rdrr.io/pkg/detectRUNS/man/plot_Runs.html));

- `stacked`: the runs of each group stacked along each chromosome, one
  run per line, ordered by start (as
  [`detectRUNS::plot_StackedRuns()`](https://rdrr.io/pkg/detectRUNS/man/plot_StackedRuns.html));

- `coverage`: the percentage of individuals of each group in a run along
  each chromosome (as
  [`detectRUNS::plot_SnpsInRuns()`](https://rdrr.io/pkg/detectRUNS/man/plot_SnpsInRuns.html)
  and
  [`detectRUNS::plot_manhattanRuns()`](https://rdrr.io/pkg/detectRUNS/man/plot_manhattanRuns.html),
  computed from the runs rather than at each locus);

- `froh`: the inbreeding coefficient F_ROH, the total length of the runs
  of an individual over the length of the genome (with
  `froh_by = "individual"`, one point per individual, by group) or of
  each chromosome (with `froh_by = "chromosome"`, the mean of each
  group) (as
  [`detectRUNS::plot_InbreedingChr()`](https://rdrr.io/pkg/detectRUNS/man/plot_InbreedingChr.html));
  the length of a chromosome is its largest position, as in
  `detectRUNS`;

- `violin`: the distribution of the total (or mean, see `method`) length
  of the runs of each individual, by group (as
  [`detectRUNS::plot_ViolinRuns()`](https://rdrr.io/pkg/detectRUNS/man/plot_ViolinRuns.html));

- `length_class`: the percentage of runs in each class of length, by
  group (as
  [`detectRUNS::plot_DistributionRuns()`](https://rdrr.io/pkg/detectRUNS/man/plot_DistributionRuns.html)).

## Usage

``` r
# S3 method for class 'indiv_roh'
autoplot(
  object,
  type = c("runs", "stacked", "coverage", "froh", "violin", "length_class"),
  chromosome = NULL,
  froh_by = c("individual", "chromosome"),
  method = c("sum", "mean"),
  classes = c(0, 2, 4, 8, 16),
  ...
)
```

## Arguments

- object:

  an object of class `indiv_roh`, as returned by
  [`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_indiv_roh.md)

- type:

  the type of plot (`runs`, `stacked`, `coverage`, `froh`, `violin` or
  `length_class`)

- chromosome:

  the chromosomes to plot (by default, all those with runs; for `froh`
  with `froh_by = "individual"`, the genome is always used)

- froh_by:

  for `froh`, whether to plot F_ROH of each individual over the genome
  (`individual`), or the mean F_ROH of each group on each chromosome
  (`chromosome`)

- method:

  for `violin`, whether to plot the total (`sum`) or the mean (`mean`)
  length of the runs of each individual

- classes:

  for `length_class`, the lower bounds of the classes of length of the
  runs, in Mb (the last class has no upper bound)

- ...:

  not currently used.

## Value

a `ggplot2` object

## Details

Individuals with no runs count in `coverage`, `froh` and `violin` (with
no runs, or a length of 0).

`autoplot` produces simple plots to quickly inspect an object. They are
not customisable; we recommend that you use `ggplot2` to produce
publication ready plots.

## Examples

``` r
# a simulated dataset of three populations with increasing inbreeding
bed_file <- system.file("extdata", "roh", "roh_sim.bed",
  package = "tidypopgen"
)
sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
sim_roh <- windows_indiv_roh(sim_gt %>% group_by(population))

autoplot(sim_roh, type = "runs", chromosome = "3")

autoplot(sim_roh, type = "coverage")

autoplot(sim_roh, type = "froh")

autoplot(sim_roh, type = "length_class")
```
