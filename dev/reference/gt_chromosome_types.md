# Types of chromosomes: autosomes, sex chromosomes and mtDNA

The loci of a `gen_tibble` are on autosomes, on sex chromosomes (X and Y
in XY systems, Z and W in ZW systems), on their pseudo-autosomal regions
(PAR), or on the mitochondrial genome (MT, which can also be used for
the genomes of other organelles, e.g. chloroplasts).
`gt_set_chromosome_types()` declares the type of chromosomes that can
not be recognised from their names (e.g. the numeric codes of PLINK),
and the pseudo-autosomal regions of the X chromosome.
`gt_chromosome_types()` lists the chromosomes with their type, and
`loci_chromosome_types()` gives the type of each locus (e.g. to select
the autosomal loci with
[`select_loci_if()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/select_loci_if.md)).

## Usage

``` r
gt_set_chromosome_types(.x, types = NULL, chr_set = NULL, par = NULL)

gt_chromosome_types(.x)

loci_chromosome_types(.x, .col = "genotypes", ...)

# S3 method for class 'tbl_df'
loci_chromosome_types(.x, .col = "genotypes", ...)

# S3 method for class 'vctrs_bigSNP'
loci_chromosome_types(.x, .col = "genotypes", ...)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  (for `loci_chromosome_types()`, also its `genotypes` column).

- types:

  a named character vector giving the type of chromosomes: its names are
  the chromosomes (as in
  [`loci_chromosomes()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_chromosomes.md);
  for chromosomes stored as integer codes, their names, see
  [`gt_lookup_chromosomes()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)),
  and its values their types, one of "autosome", "X", "Y", "Z", "W",
  "PAR" and "MT". For example, `c("23" = "X", "24" = "Y")`.

- chr_set:

  the number of autosomes, as in PLINK's `--chr-set`, to declare the
  numeric codes of X, Y, XY and MT (see details), or "human". The types
  given in `types` take precedence.

- par:

  the pseudo-autosomal regions of the X chromosome: either a build of
  the human genome ("b36", "b37", "b38", or their names "hg18", "hg19",
  "hg38"), whose boundaries apply to the chromosome of type "X", or a
  data frame with columns `chromosome`, `start` and `end` (positions,
  both included) giving the regions of each chromosome of type "X" that
  are pseudo-autosomal.

- .col:

  the column to be used when a tibble (or grouped tibble is passed
  directly to the function). This defaults to "genotypes" and can only
  take that value. There is no need for the user to set it, but it is
  included to resolve certain tidyselect operations.

- ...:

  currently unused.

## Value

`gt_set_chromosome_types()` returns the `gen_tibble` with the declared
types; `gt_chromosome_types()` returns a tibble with columns
`chromosome`, `type` and `n_loci` (the number of loci), with one row for
each chromosome and type (an X chromosome with pseudo-autosomal regions
has two rows); `loci_chromosome_types()` returns a factor with the type
of each locus.

## Details

Chromosomes are recognised from their names, regardless of case and of a
prefix "chr", "chr\_", "chromosome" or "chromosome\_" (e.g. "X", "chrX"
or "chromosome_x"): "X", "Y", "Z" and "W" are sex chromosomes; "XY",
"PAR", "PAR1" and "PAR2" are pseudo-autosomal regions; and "MT", "M" and
"0M" are the mitochondrial genome. All other chromosomes are autosomes.
Numbers are never recognised, as their meaning depends on the species
(PLINK codes the human X as 23, but chromosome 23 is an autosome in a
dog): declare them with `types` or `chr_set`.

The types declared with `gt_set_chromosome_types()` are kept with the
`gen_tibble` (they are saved by
[`gt_save()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_save.md));
each call replaces the previous declaration, and a call without `types`,
`chr_set` and `par` removes it. A declaration can also bring a
recognised chromosome back among the autosomes (e.g.
`types = c(X = "autosome")`).

Functions that pool loci and assume that they are on diploid autosomes
only use the autosomal loci and the pseudo-autosomal regions, as PLINK
does (e.g. for KING), and say how many loci they left out: the
individual statistics
[`indiv_het_obs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/indiv_het_obs.md),
[`indiv_inbreeding()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/indiv_inbreeding.md)
and
[`qc_report_indiv()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/qc_report_indiv.md);
the pairwise statistics
[`pairwise_king()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_king.md),
[`pairwise_ibs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_ibs.md),
[`pairwise_dist()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_dist.md),
[`pairwise_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_allele_sharing.md)
and
[`pairwise_grm()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_grm.md);
[`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md)
(and the deprecated `gt_pca_*SVD()`) and
[`gt_pcadapt()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pcadapt.md);
[`gt_extract_f2()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_extract_f2.md),
[`gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_admixture.md)
and
[`gt_snmf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_snmf.md);
[`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_indiv_roh.md);
and the population statistics
([`pop_fst()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_fst.md),
[`pop_fis()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_fis.md),
[`pop_het_obs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_het_obs.md),
[`pop_het_exp()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_het_exp.md),
[`pop_gene_div()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_het_exp.md),
[`pop_global_stats()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_global_stats.md),
[`pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_pop_fst.md),
[`nwise_pop_pbs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/nwise_pop_pbs.md),
[`pop_tajimas_d()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_tajimas_d.md),
the site frequency spectra
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_sfs.md),
[`pairwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_pop_sfs.md),
[`nwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/nwise_pop_sfs.md),
[`bootstrap_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/bootstrap_sfs.md),
and their `windows_*()` and `jackknife_*()` versions), whose results by
locus or by window then cover the autosomal loci only.
[`gt_pseudohaploid()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pseudohaploid.md)
also tests autosomal loci. The functions that return one value per locus
(`loci_*()`,
[`qc_report_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/qc_report_loci.md)),
LD clumping and pruning, and imputation use all loci, without treating
the sex chromosomes of males as haploid. To include a sex chromosome or
the mitochondrial genome in the analyses that use autosomes only,
declare it as an autosome.

PLINK codes chromosomes with numbers: with n autosomes, X is n+1, Y is
n+2, the pseudo-autosomal region of X ("XY") is n+3, and MT is n+4 (23
to 26 in humans, as set by `--chr-set n`). `chr_set = n` declares those
codes; `chr_set = "human"` is `chr_set = 22`, plus the codes of
EIGENSOFT for the mitochondrial genome (90) and XY (91).

PLINK's `--split-x` moves the loci of the pseudo-autosomal regions of
the human X chromosome to the chromosome "XY", using the boundaries of a
build of the human genome: a locus on X is in a pseudo-autosomal region
if its position is at most the last position of PAR1, or at least the
first position of PAR2. `par` gives the type "PAR" to those loci without
changing their chromosome: "b36" (or "hg18") uses 2709521 and 154584237,
"b37" (or "hg19") 2699520 and 154931044, and "b38" (or "hg38") 2781479
and 155701383, as PLINK.

## References

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D,
Maller J, Sklar P, de Bakker PIW, Daly MJ, Sham PC (2007) PLINK: a tool
set for whole-genome association and population-based linkage analyses.
American Journal of Human Genetics, 81(3), 559-575.
[doi:10.1086/519795](https://doi.org/10.1086/519795)

## Examples

``` r
bed_file <- system.file("extdata", "pop_a.bed", package = "tidypopgen")
example_gt <- gen_tibble(bed_file,
  backingfile = tempfile("pop_a_"),
  quiet = TRUE
)

# chromosome 23 is taken as an autosome
gt_chromosome_types(example_gt)
#> # A tibble: 4 × 3
#>   chromosome type     n_loci
#>   <chr>      <fct>     <int>
#> 1 1          autosome     11
#> 2 2          autosome      3
#> 3 3          autosome      1
#> 4 23         autosome      1

# declare the human codes of PLINK: 23 is the X chromosome
example_gt <- gt_set_chromosome_types(example_gt, chr_set = "human")
gt_chromosome_types(example_gt)
#> # A tibble: 4 × 3
#>   chromosome type     n_loci
#>   <chr>      <fct>     <int>
#> 1 1          autosome     11
#> 2 2          autosome      3
#> 3 3          autosome      1
#> 4 23         X             1
loci_chromosome_types(example_gt)
#>  [1] autosome autosome autosome autosome autosome autosome autosome autosome
#>  [9] autosome autosome autosome autosome autosome autosome autosome X       
#> Levels: autosome X Y Z W PAR MT

# keep the autosomal loci only
example_gt %>%
  select_loci_if(loci_chromosome_types(genotypes) == "autosome") %>%
  count_loci()
#> [1] 15
```
