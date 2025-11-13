# PLINK cheatsheet

The following flags and corresponding tidypopgen functions are based on
plink version 1.9.

## File management and reading data:

| PLINK 1.9                   | `tidypopgen`                                                                                                           |
|-----------------------------|------------------------------------------------------------------------------------------------------------------------|
| –make-bed –out              | `gt_as_plink(data, file = my_file, type = "bed")`                                                                      |
| –recode                     | `gt_as_plink(data, file = my_file, type = "ped")`                                                                      |
| –recode vcf                 | `gt_as_vcf(data, file = my_file)`                                                                                      |
| –allele1234 and –alleleACGT | See [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md) parameter ‘valid_alleles’ |

PLINK flags –update-alleles, –allele1234, and –alleleACGT, all alter the
coding of alleles. In `tidypopgen`, valid alleles are supplied when
reading in a `gen_tibble`.

## Quality control:

| PLINK            | `tidypopgen`                   |
|------------------|--------------------------------|
| –maf             | `data %>% loci_maf()`          |
| –geno            | `data %>% loci_missingness()`  |
| –hwe             | `data %>% loci_hwe()`          |
| –freq –a2-allele | `data %>% loci_alt_freq()`     |
| –mind            | `data %>% indiv_missingness()` |
| –het             | `data %>% indiv_het_obs()`     |

To filter out variants in `tidypopgen`, in a similar way to PLINK flags
such as –extract or –autosome, it is necessary to use the `gen_tibble`
with
[`select_loci_if()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/select_loci_if.md).
For example:

``` r
data %>% select_loci_if(loci_chromosomes(genotypes) %in% c(1:22))
#> # A gen_tibble: 15 loci
#> # A tibble:     5 × 4
#>   id    population sex    genotypes
#>   <chr> <chr>      <fct> <vctr_SNP>
#> 1 GRC24 pop_a      male   [0,0,...]
#> 2 GRC25 pop_a      male   [0,0,...]
#> 3 GRC26 pop_a      male   [0,1,...]
#> 4 GRC27 pop_a      male   [0,0,...]
#> 5 GRC28 pop_a      male   [0,0,...]
```

will select autosomal loci in the same way as –autosome. Or
alternatively:

``` r
my_snps <- c("rs4477212", "rs3094315", "rs3131972", "rs12124819", "rs11240777")

data %>%
  select_loci_if(loci_names(genotypes) %in% my_snps) %>%
  show_loci()
#> # A tibble: 5 × 7
#>   big_index name       chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr>      <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         1 rs4477212  1             82154            0 A          NA        
#> 2         2 rs3094315  1            752566            0 A          G         
#> 3         3 rs3131972  1            752721            0 G          A         
#> 4         4 rs12124819 1            776546            0 A          NA        
#> 5         5 rs11240777 1            798959            0 G          A
```

will select loci from a previously defined set in the same way as
–extract.

Similarly, to filter out individuals, as might be performed with –keep
in PLINK, requires using filter:

``` r
my_individuals <- c("GRC14300079", "GRC14300142", "GRC14300159")

data %>% filter(id %in% my_individuals)
#> # A gen_tibble: 16 loci
#> # A tibble:     0 × 4
#> # ℹ 4 variables: id <chr>, population <chr>, sex <fct>, genotypes <vctr_SNP>
```

## Handling linkage

Linkage disequilibrium is managed through clumping in `tidypopgen` with
[`loci_ld_clump()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_ld_clump.md).

This option is similar to the –indep-pairwise flag in PLINK, but results
in a more even distribution of loci when compared to LD pruning.

To explore why clumping is preferable to pruning, see
<https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html>

## Quality control for relatedness (KING)

| KING       | `tidypopgen`                                                                                                       |
|------------|--------------------------------------------------------------------------------------------------------------------|
| –kinship   | [`pairwise_king()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_king.md)                     |
| –distance  | [`pairwise_ibs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_ibs.md)                       |
| –unrelated | [`filter_high_relatedness()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/filter_high_relatedness.md) |

[`pairwise_king()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_king.md)
implements the KING-robust estimator of kinship, equivalent to –kinship
in KING. To remove related individuals, the user can pass the kinship
matrix of
[`pairwise_king()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_king.md)
and a relatedness threshold (a numeric KING kinship coefficient) to
[`filter_high_relatedness()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/filter_high_relatedness.md),
which will return the largest possible set of individuals with no
relationships above the threshold.

[`pairwise_king()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_king.md)
also forms part of
[`qc_report_indiv()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/qc_report_indiv.md)
through the parameter `kings_threshold`. To remove related individuals
in
[`qc_report_indiv()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/qc_report_indiv.md),
the user can pass either a relatedness threshold, or a string of either
“first” or “second” to remove any first degree or second degree
relationships from the dataset. This second option is similar to using
–unrelated –degree 1 or –unrelated –degree 2 in KING.
[`qc_report_indiv()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/qc_report_indiv.md)
returns a dataframe including columns ‘id’ and ‘to_keep’, showing the ID
of individuals and a logical column of whether the individual should be
removed to retain the largest possible set of individuals with no
relationships above the threshold.

## Merging datasets:

| PLINK      | `tidypopgen`                                                                                   |
|------------|------------------------------------------------------------------------------------------------|
| –bmerge    | [`rbind()`](https://rdrr.io/r/base/cbind.html)                                                 |
| –flip-scan | [`rbind_dry_run()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/rbind_dry_run.md) |
| –flip      | use ‘flip_strand = TRUE’ in [`rbind()`](https://rdrr.io/r/base/cbind.html)                     |

In PLINK, data merging can fail due to strand inconsistencies that are
not addressed prior to merging. PLINK documentation suggests to users to
try a ‘trial flip’ of data to address this, and then to ‘unflip’ any
errors that remain. In `tidypopgen`, when data are merged with rbind,
strand inconsistencies are identified and automatically flipped,
avoiding multiple rounds of flipping before merging.

PLINK does allow users to identify inconsistencies prior to merging with
–flip-scan, and this functionality is included in the `tidypopgen`
rbind_dry_run(). rbind_dry_run() reports the numeric overlap of
datasets, alongside the number of SNPs to ‘flip’ in the new target
dataset, as well as the number of ambiguous SNPs.

Data are only merged one set at a time, there is no equivalent to
–merge-list.

## Analysis:

| PLINK    | `tidypopgen`                                                                                                                                                                  |
|----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| –pca     | See [`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md) for pca options                                                                          |
| –fst     | [`pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_pop_fst.md) with [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) |
| –homozyg | [`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_indiv_roh.md)                                                                        |
