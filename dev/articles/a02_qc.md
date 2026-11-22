# Quality Control

## Quality control for SNP datasets

tidypopgen has two key functions to examine the quality of data across
loci and across individuals: `qc_report_loci` and `qc_report_indiv`.
This vignette uses a simulated data set to illustrate these methods of
data cleaning.

If your `gen_tibble` includes several populations, group it by
population with
[`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) as
soon as you have read it, before any quality control: `qc_report_loci`
then tests Hardy-Weinberg equilibrium within each population, and
`qc_report_indiv` estimates relatedness within each population (see
[Grouping data](#grouping-data) below). The simulated data set that we
use here is a single population.

## Read data into gen_tibble format

``` r

library(tidypopgen)
```

    ## Loading required package: dplyr

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Loading required package: tibble

``` r

data <- gen_tibble(
  system.file("extdata/related/families.bed",
    package = "tidypopgen"
  ),
  quiet = TRUE, backingfile = tempfile(),
  valid_alleles = c("1", "2")
)
```

First, lets take a look at the data we have.

``` r

count_loci(data)
```

    ## [1] 961

``` r

nrow(data)
```

    ## [1] 12

We can see that there are 961 SNPs and 12 individuals in this dataset.

## Quality control for individuals

``` r

individual_report <- qc_report_indiv(data)
summary(individual_report)
```

    ##     het_obs        missingness     
    ##  Min.   :0.3467   Min.   :0.03018  
    ##  1st Qu.:0.3736   1st Qu.:0.03642  
    ##  Median :0.3790   Median :0.03798  
    ##  Mean   :0.3799   Mean   :0.03867  
    ##  3rd Qu.:0.3901   3rd Qu.:0.04110  
    ##  Max.   :0.4015   Max.   :0.04579

The output of `qc_report_indiv` supplies observed heterozygosity per
individual, and rate of missingness per individual as standard.

These data can also be visualised using autoplot:

``` r

autoplot(individual_report)
```

![Scatter plot of missingness proportion and observed heterozygosity for
each individual](a02_qc_files/figure-html/unnamed-chunk-4-1.png)

Here, the red line indicates a threshold for proportion of missing loci,
which is set as 5% by default, and can be altered using the
`miss_threshold` argument. Similarly, the blue lines are drawn at 2
(inner) and 3 (outer) standard deviations from the mean observed
heterozygosity. These thresholds can be used to visualise outliers and
consider how to filter datasets.

We can see above that most individuals have low missingness, none are
above the default 5% threshold. However, if we wanted to filter
individuals to remove those with more than 4.5% of their genotypes
missing, we can use `filter`. The report has one row per individual, in
the same order as `data`, so we take the ids of the individuals to keep
and filter on them; this also works on a grouped `gen_tibble`, where
`filter` works on one group at a time.

``` r

keep_ids <- data$id[individual_report$missingness < 0.045]
data <- data %>% filter(id %in% keep_ids)
nrow(data)
```

    ## [1] 10

And if we wanted to remove outliers with particularly high or low
heterozygosity, we can again do so by using `filter`. As an example,
here we remove observations that lie more than 2 standard deviations
from the mean.

``` r

mean_val <- mean(individual_report$het_obs)
sd_val <- stats::sd(individual_report$het_obs)

lower <- mean_val - 2 * (sd_val)
upper <- mean_val + 2 * (sd_val)

data <- data %>% filter(indiv_het_obs(genotypes) > lower)
data <- data %>% filter(indiv_het_obs(genotypes) < upper)
nrow(data)
```

    ## [1] 9

Next, we can look at relatedness within our sample. If the parameter
`kings_threshold` is provided to
[`qc_report_indiv()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/qc_report_indiv.md),
then the report also calculates a KING coefficient of relatedness matrix
using the sample. The `kings_threshold` is used to provide an output of
the largest possible group with no related individuals in the third
column `to_keep`. This boolean column recommends which individuals to
remove (FALSE) and to keep (TRUE) to achieve an unrelated sample.

For example, we can use `kings_threshold = "first"` to specify that we
want to remove first degree relatives.

``` r

individual_report <- qc_report_indiv(data, kings_threshold = "first")
summary(individual_report)
```

    ##     het_obs        missingness       to_keep                id   
    ##  Min.   :0.3688   Min.   :0.03018   Mode :logical   Length   :9  
    ##  1st Qu.:0.3774   1st Qu.:0.03642   FALSE:1         N.unique :9  
    ##  Median :0.3795   Median :0.03746   TRUE :8         N.blank  :0  
    ##  Mean   :0.3851   Mean   :0.03735                   Min.nchar:1  
    ##  3rd Qu.:0.3985   3rd Qu.:0.04058                   Max.nchar:2  
    ##  Max.   :0.4015   Max.   :0.04266

We can then remove the recommended individuals by using:

``` r

data <- data %>%
  filter(id %in% individual_report$id[individual_report$to_keep])
```

We can now view a summary of our cleaned data set again, showing that
our data has reduced from the original; 12 individuals, to 8.

``` r

summary(data)
```

    ##          id        genotypes 
    ##  Length   :8   Length   : 8  
    ##  N.unique :8   N.unique : 1  
    ##  N.blank  :0   N.blank  : 0  
    ##  Min.nchar:1   Min.nchar:16  
    ##  Max.nchar:2   Max.nchar:16

## Quality control for loci

``` r

loci_report <- qc_report_loci(data)
```

    ## This gen_tibble is not grouped. For Hardy-Weinberg equilibrium, `qc_report_loci()` will assume individuals are part of the same population and HWE test p-values will be calculated across all individuals. If you wish to calculate HWE p-values within populations or groups, please use`group_by()` before calling `qc_report_loci()`.

``` r

summary(loci_report)
```

    ##        snp_id         maf          missingness         hwe_p        
    ##  Length   :961   Min.   :0.0000   Min.   :0.0000   Min.   :0.01538  
    ##  N.unique :961   1st Qu.:0.1429   1st Qu.:0.0000   1st Qu.:0.32867  
    ##  N.blank  :  0   Median :0.2500   Median :0.0000   Median :0.53333  
    ##  Min.nchar:  1   Mean   :0.2648   Mean   :0.0372   Mean   :0.50413  
    ##  Max.nchar:  3   3rd Qu.:0.3750   3rd Qu.:0.1250   3rd Qu.:0.69231  
    ##                  Max.   :0.5000   Max.   :0.3750   Max.   :0.76503

The output of `qc_report_loci` supplies minor allele frequency, rate of
missingness, and a Hardy-Weinberg exact p-value for each SNP. These data
can be visualised in autoplot :

``` r

autoplot(loci_report, type = "overview")
```

![UpSet plot giving counts of snps over the threshold for: missingness,
minor allele frequency, and Hardy-Weinberg equilibrium
P-value](a02_qc_files/figure-html/unnamed-chunk-11-1.png)

Using ‘overview’ provides an Upset plot. Upset plots are designed to
show the intersection of different sets in the same way as a Venn
diagram. SNPs can be divided into ‘sets’ that each pass predefined
quality control threshold; a set of SNPs with missingness under a given
threshold, a set of SNPs with MAF above a given threshold, and a set of
SNPs with a Hardy-Weinberg exact p-value that falls above a given
significance level.

The thresholds for each parameter, (percentage of missingness that is
accepted, minor allele frequency cutoff, and Hardy-Weinberg equilibrium
p-value) can be adjusted using the parameters provided in autoplot. For
example, lets adjust the thresholds to 5% missingness, 10% minor allele
frequency, and a Hardy-Weinberg p-value of 0.05.

``` r

autoplot(loci_report,
  type = "overview",
  miss_threshold = 0.05,
  maf_threshold = 0.1,
  hwe_p = 0.05
)
```

![Upset plot as above, with adjusted
thresholds](a02_qc_files/figure-html/unnamed-chunk-12-1.png)

The upset plot then visualises our 961 SNPs within their respective
sets. The number above the first bar indicates that 607 SNPs occur in
all 3 sets, meaning 607 SNPs pass all of our QC thresholds. The combined
total of the first and second bars represents the number of SNPs that
pass our MAF and HWE thresholds, here 804 SNPs.

To examine each QC measure in further detail, we can plot a different
summary panel.

``` r

autoplot(loci_report,
  type = "all",
  miss_threshold = 0.03,
  maf_threshold = 0.02,
  hwe_p = 0.01
)
```

![Four panel plot, containing: a histogram of the proportion of missing
data for snps with minor allele frequency above the threshold, a
histogram of the proportion of missing data for snps with minor allele
freqency below the threshold, a histogram of HWE exact test p-values,
and a histogram of significant HWE exact test
p-values](a02_qc_files/figure-html/unnamed-chunk-13-1.png)

We can then begin to consider how to quality control this raw data set.
The loci report already contains the minor allele frequency (`maf`), the
rate of missingness (`missingness`) and the Hardy-Weinberg exact p-value
(`hwe_p`) of every SNP, so we can use it to decide which SNPs to keep
without calculating these statistics again. First, let’s look at each
measure in turn to choose our thresholds.

Let’s start with the minor allele frequency. We can visualise the MAF
distribution using:

``` r

autoplot(loci_report, type = "maf")
```

![Histogram of minor allele
frequency](a02_qc_files/figure-html/unnamed-chunk-14-1.png)

Here we can see there are some monomorphic SNPs in the data set. We will
remove loci with a minor allele frequency lower than 2%, which is
equivalent to using plink –maf 0.02.

Next, we can look at the rate of missingness. Lets say we want to remove
SNPs that are missing in more than 5% of individuals, equivalent to
using plink –geno 0.05

``` r

autoplot(loci_report, type = "missing", miss_threshold = 0.05)
```

![Histogram of the proportion of missing
data](a02_qc_files/figure-html/unnamed-chunk-15-1.png)

We can see here that most SNPs have low missingness, under our 5%
threshold, some do, however, have missingness over our threshold.

Finally, we may want to remove SNPs that show significant deviation from
Hardy-Weinberg equilibrium, if our study design requires. To visualise
SNPs with significant p-values in the Hardy-Weinberg exact test, we can
again call autoplot:

``` r

autoplot(loci_report, type = "significant hwe", hwe_p = 0.01)
```

![Histogram of significant HWE exact test
p-values](a02_qc_files/figure-html/unnamed-chunk-16-1.png)

None of the SNPs in our data are significant, however there may be
circumstances where we would want to cut out the most extreme cases, if
these data were real, these cases could indicate genotyping errors.

Now that we have chosen our thresholds, we can apply all three in one
go, using the `loci_report` we have just created. Each column of the
report has one value per SNP, in the same order as the loci in our
`gen_tibble`. We can therefore combine the three conditions into a
single logical vector, which is `TRUE` for the SNPs that we want to
keep, and pass it to `select_loci_if`:

``` r

loci_to_keep <- loci_report$maf > 0.02 &
  loci_report$missingness < 0.05 &
  loci_report$hwe_p > 0.01

data <- data %>% select_loci_if(loci_to_keep)
count_loci(data)
```

    ## [1] 693

Note that we need to apply all of the filters at once. The `loci_report`
was calculated on a `gen_tibble` with 961 SNPs, so it has one row for
each of those 961 SNPs. As soon as `select_loci_if` has removed some
SNPs, `data` has fewer loci than the `loci_report` has rows, and the two
no longer match: a logical vector built from the old report
(e.g. `loci_report$missingness < 0.05`) would have the wrong length for
the reduced `gen_tibble`. For the same reason, a report cannot be reused
to filter a `gen_tibble` that has since been reduced. The report also
describes the data as it was when it was created, so its values would be
out of date for the reduced data anyway.

If we want to apply the filters in separate steps, there are two
options. We can calculate the statistics again within `select_loci_if`,
which always operates on the current loci of the `gen_tibble`:

``` r

data <- data %>% select_loci_if(loci_maf(genotypes) > 0.02)
data <- data %>% select_loci_if(loci_missingness(genotypes) < 0.05)
data <- data %>% select_loci_if(loci_hwe(genotypes) > 0.01)
```

Or we can create a new report with
[`qc_report_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/qc_report_loci.md)
after each step, so that it matches the reduced `gen_tibble`. On a
grouped `gen_tibble`, only the second option tests Hardy-Weinberg
equilibrium within each population: `select_loci_if` evaluates its
condition on all the individuals, ignoring the groups, so
`loci_hwe(genotypes)` within it tests all the individuals together (see
[Grouping data](#grouping-data)).

## Linkage Disequilibrium

For further analyses, it may be necessary to control for linkage in the
data set. tidypopgen provides LD clumping, with
[`loci_ld_clump()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_ld_clump.md).
This option is similar to the –indep-pairwise flag in plink, but results
in a more even distribution of loci when compared to LD pruning.
Pruning, as plink –indep-pairwise, is available with
[`loci_ld_prune()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_ld_prune.md).

To explore why clumping is preferable to pruning, see
<https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html>

Missing genotypes do not need to be imputed: the correlation between two
loci is computed over the individuals genotyped at both (as in PLINK).

In this example, if we want to remove SNPs with a correlation greater
than 0.2 within windows of 10 kb, we can set these parameters with
`thr_r2` and `size` respectively (with `use_positions = FALSE`, `size`
would be a number of SNPs instead).

``` r

to_keep_ld <- loci_ld_clump(data, thr_r2 = 0.2, size = 10)
head(to_keep_ld)
```

    ## [1] FALSE FALSE FALSE FALSE  TRUE FALSE

`loci_ld_clump` provides a boolean vector the same length as our list of
SNPs, telling us which to keep in the data set. We can then use this
list to create a pruned version of our data:

``` r

ld_data <- data %>%
  select_loci_if(loci_ld_clump(genotypes, thr_r2 = 0.2, size = 10))
```

Pruning as `plink --indep-pairwise 10 1 0.2` (a window of 10 SNPs, moved
by 1 SNP at a time) keeps this number of SNPs:

``` r

data %>%
  loci_ld_prune(window_size = 10, step_size = 1, thr_r2 = 0.2) %>%
  sum()
```

    ## [1] 203

## Save

The benefit of operating on a `gen_tibble` is that each quality control
step can be observed visually, and easily reversed if necessary.

When we are happy with the quality of our data, we can create and save a
final quality controlled version of our `gen_tibble` using `gt_save`.

``` r

gt_save(ld_data, file_name = tempfile())
```

    ## 
    ## gen_tibble saved to /tmp/RtmpJlRgOk/file2a5d9e8b711.gt

    ## using FBM RDS: /tmp/RtmpJlRgOk/file2a5d1076691c.rds

    ## with FBM backing file: /tmp/RtmpJlRgOk/file2a5d1076691c.bk

    ## make sure that you do NOT delete those files!

    ## to reload the gen_tibble in another session, use:

    ## gt_load('/tmp/RtmpJlRgOk/file2a5d9e8b711.gt')

    ## [1] "/tmp/RtmpJlRgOk/file2a5d9e8b711.gt"  
    ## [2] "/tmp/RtmpJlRgOk/file2a5d1076691c.rds"
    ## [3] "/tmp/RtmpJlRgOk/file2a5d1076691c.bk"

## Grouping data

For some quality control measures, if you have a gen_tibble that
includes multiple datasets you may want to group by population before
running the quality control. This can be done using `group_by`.

First, lets add some imaginary population data to our gen_tibble:

``` r

data <- data %>% mutate(population = c(rep("A", 4), rep("B", 4)))
```

We can then group by population:

``` r

grouped_data <- data %>% group_by(population)
```

and run quality control on each group:

``` r

grouped_loci_report <- grouped_data %>%
  qc_report_loci()
head(grouped_loci_report)
```

    ## # A tibble: 6 × 4
    ##   snp_id    maf missingness hwe_p
    ##   <chr>   <dbl>       <dbl> <dbl>
    ## 1 2      0.375            0 0.429
    ## 2 3      0.188            0 1    
    ## 3 5      0.312            0 1    
    ## 4 8      0.25             0 1    
    ## 5 10     0.5              0 1    
    ## 6 11     0.0625           0 1

The loci report that we receive here will calculate Hardy-Weinberg
equilibrium for each SNP within each population separately, providing a
Bonferroni corrected p-value for each SNP.

Similarly, we can run a quality control report for individuals within
each population:

``` r

grouped_individual_report <- grouped_data %>%
  qc_report_indiv(kings_threshold = 0.177)
head(grouped_individual_report)
```

    ## # A tibble: 6 × 5
    ##   het_obs missingness id    group to_keep
    ##     <dbl>       <dbl> <chr> <chr> <lgl>  
    ## 1   0.380           0 1     A     TRUE   
    ## 2   0.375           0 4     A     TRUE   
    ## 3   0.390           0 5     A     TRUE   
    ## 4   0.380           0 6     A     TRUE   
    ## 5   0.410           0 7     B     TRUE   
    ## 6   0.378           0 8     B     TRUE

This is important when we have related individuals, as background
population structure can affect the filtering of relatives.

## Grouped functions

It is also possible to run `loci` and `indiv` functions on grouped data.
This is useful when you want to run the same quality control on each
group of data, but don’t want to split the data into separate
gen_tibbles.

``` r

loci_maf_grouped <- grouped_data %>%
  loci_maf()
```

Grouped functions are built for efficiency and surpass the use of
applying a function with `group_map`.

Note that `select_loci_if` ignores the groups: all populations share the
same loci, so a locus is kept or removed for all of them, and
`select_loci_if(loci_maf(genotypes) > 0.02)` uses the minor allele
frequency over all individuals, even on a grouped `gen_tibble` (as
`qc_report_loci` does for the minor allele frequency and missingness).
To filter on the frequencies within populations, we combine the results
of the grouped `loci_maf` into one value per locus. For example, to keep
the loci with a minor allele frequency above 2% in every population:

``` r

maf_by_pop <- grouped_data %>% loci_maf(type = "matrix")
grouped_data_maf <- grouped_data %>%
  select_loci_if(apply(maf_by_pop, 1, min) > 0.02)
count_loci(grouped_data_maf)
```

    ## [1] 534
