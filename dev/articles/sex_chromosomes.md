# Sex chromosomes and mtDNA

Most analyses in population genetics assume diploid loci, but the ploidy
of a locus on a sex chromosome depends on the sex of the individual. In
XY systems (e.g. mammals), males have one X and one Y, and females two X
and no Y; in ZW systems (e.g. birds), females have one Z and one W. The
mitochondrial genome (mtDNA) is haploid in everyone. Only the
pseudo-autosomal regions (PAR) of the sex chromosomes behave like
autosomes. This article shows how `tidypopgen` deals with these loci:

1.  how it recognises sex chromosomes and mtDNA;
2.  how it stores the genotypes of hemizygous individuals;
3.  how to check the recorded sex of individuals;
4.  which analyses use autosomes only;
5.  how sex chromosomes are written to PLINK and VCF files.

`tidypopgen` follows PLINK’s conventions throughout.

## The data

We simulate a small dataset of 10 females and 10 males, with 50 loci on
an autosome, 100 on the X chromosome, 10 on Y and 5 on mtDNA. As in
PLINK, the genotypes of males on X and Y (and of everyone on mtDNA),
which have a single allele, are stored as homozygotes (0 or 2).

``` r

library(tidypopgen)
set.seed(42)
n_females <- 10
n_males <- 10
male <- rep(c(FALSE, TRUE), c(n_females, n_males))
# diploid genotypes, and haploid ones stored as homozygotes
diploid <- function(p, n) t(replicate(n, rbinom(length(p), 2, p)))
haploid <- function(p, n) 2 * t(replicate(n, rbinom(length(p), 1, p)))
p_x <- runif(100, 0.1, 0.9)
x <- matrix(NA, length(male), 100)
x[!male, ] <- diploid(p_x, n_females)
x[male, ] <- haploid(p_x, n_males)
y <- matrix(NA, length(male), 10)
y[male, ] <- haploid(runif(10, 0.2, 0.8), n_males)
geno <- cbind(
  diploid(runif(50, 0.1, 0.9), length(male)),
  x,
  y,
  haploid(runif(5, 0.2, 0.8), length(male))
)
loci <- data.frame(
  name = paste0("snp", 1:165),
  chromosome = rep(c("1", "X", "Y", "MT"), c(50, 100, 10, 5)),
  position = c(1:50, 1:100, 1:10, 1:5) * 1000,
  genetic_dist = 0,
  allele_ref = "A",
  allele_alt = "G"
)
sex_gt <- gen_tibble(geno,
  loci = loci,
  indiv_meta = data.frame(
    id = paste0("ind", 1:20),
    population = "pop1",
    sex = ifelse(male, "male", "female")
  ),
  backingfile = tempfile("sex_chromosomes_"),
  quiet = TRUE
)
sex_gt
#> # A gen_tibble: 165 loci
#> # A tibble:     20 × 4
#>    id    population sex     genotypes
#>    <chr> <chr>      <chr>  <vctr_SNP>
#>  1 ind1  pop1       female  [1,2,...]
#>  2 ind2  pop1       female  [1,2,...]
#>  3 ind3  pop1       female  [1,1,...]
#>  4 ind4  pop1       female  [1,1,...]
#>  5 ind5  pop1       female  [2,1,...]
#>  6 ind6  pop1       female  [1,1,...]
#>  7 ind7  pop1       female  [1,0,...]
#>  8 ind8  pop1       female  [1,0,...]
#>  9 ind9  pop1       female  [0,0,...]
#> 10 ind10 pop1       female  [1,1,...]
#> 11 ind11 pop1       male    [0,1,...]
#> 12 ind12 pop1       male    [0,1,...]
#> 13 ind13 pop1       male    [2,1,...]
#> 14 ind14 pop1       male    [2,1,...]
#> 15 ind15 pop1       male    [1,1,...]
#> 16 ind16 pop1       male    [1,1,...]
#> 17 ind17 pop1       male    [2,2,...]
#> 18 ind18 pop1       male    [2,0,...]
#> 19 ind19 pop1       male    [2,1,...]
#> 20 ind20 pop1       male    [0,0,...]
```

The sex of individuals is stored in the column `sex`. It is read from
PLINK `.fam` and `.ped` files and from EIGENSOFT `.ind` files; for VCF
files, which do not record it, add it with
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) or a
join with the metadata.

## Recognising sex chromosomes and mtDNA

Chromosomes named “X”, “Y”, “Z”, “W”, “XY”, “PAR1”, “PAR2”, “MT” or “M”
(in any case, with or without “chr”) are recognised from their names:

``` r

gt_chromosome_types(sex_gt)
#> # A tibble: 4 × 3
#>   chromosome type     n_loci
#>   <chr>      <fct>     <int>
#> 1 1          autosome     50
#> 2 MT         MT            5
#> 3 X          X           100
#> 4 Y          Y            10
```

Numbers are never recognised, as their meaning depends on the species:
PLINK codes the human X chromosome as 23, but chromosome 23 is an
autosome in many other species. Declare them with
[`gt_set_chromosome_types()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_chromosome_types.md),
e.g. with PLINK’s human codes (`chr_set = "human"`, i.e. 23 to 26 for X,
Y, XY and MT), the codes of another number of autosomes (`chr_set = 38`
for a dog), or chromosome by chromosome (`types = c("23" = "X")`). The
pseudo-autosomal regions of the human X chromosome are declared with
`par = "b37"` (or “b36” and “b38”), with the boundaries of PLINK’s
`--split-x`.

[`loci_chromosome_types()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_chromosome_types.md)
gives the type of each locus, e.g. to keep the autosomes only:

``` r

sex_gt %>%
  select_loci_if(loci_chromosome_types(genotypes) == "autosome") %>%
  count_loci()
#> [1] 50
```

## Reading hemizygous genotypes

In a VCF, the genotype of a male on X is often haploid (e.g. `1` rather
than `1/1`).
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
[`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md)
and
[`gt_concat_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_concat_vcfs.md)
store a haploid genotype of a diploid individual as a homozygote (0 or
2), as PLINK does, so that the same data give the same genotypes whether
they are read from a VCF or a PLINK file. Haploid genotypes on
chromosomes recognised as sex chromosomes or mtDNA do not set the ploidy
of individuals: a VCF that starts with the X chromosome still gives
diploid males.

## Checking the sex of individuals

Males have a single X, so all their genotypes on X are homozygous, and
their inbreeding coefficient F on X is close to 1, whereas that of
females is close to 0.
[`indiv_sex_check()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/indiv_sex_check.md)
estimates the sex of each individual from its F, as PLINK’s
`--check-sex`, and compares it with the recorded sex. Let us record two
individuals with the wrong sex:

``` r

sex_gt$sex[c(1, 20)] <- c("male", "female")
indiv_sex_check(sex_gt)
#> # A tibble: 20 × 5
#>    id    sex    snp_sex status        f
#>    <chr> <fct>  <fct>   <chr>     <dbl>
#>  1 ind1  male   female  PROBLEM -0.0251
#>  2 ind2  female NA      PROBLEM  0.224 
#>  3 ind3  female female  OK       0.113 
#>  4 ind4  female female  OK       0.0580
#>  5 ind5  female female  OK      -0.0805
#>  6 ind6  female female  OK      -0.108 
#>  7 ind7  female female  OK       0.0580
#>  8 ind8  female female  OK       0.0580
#>  9 ind9  female female  OK      -0.0251
#> 10 ind10 female female  OK      -0.0251
#> 11 ind11 male   male    OK       1     
#> 12 ind12 male   male    OK       1     
#> 13 ind13 male   male    OK       1     
#> 14 ind14 male   male    OK       1     
#> 15 ind15 male   male    OK       1     
#> 16 ind16 male   male    OK       1     
#> 17 ind17 male   male    OK       1     
#> 18 ind18 male   male    OK       1     
#> 19 ind19 male   male    OK       1     
#> 20 ind20 female male    PROBLEM  1
```

The two individuals are flagged as “PROBLEM”, and so is `ind2`: with
only 100 loci on X, its F falls between the largest F of a female (0.2)
and the smallest F of a male (0.8), so its sex is unknown. With
`ycount = TRUE`, the number of genotypes on Y is also used (females
should have none).

``` r

indiv_sex_check(sex_gt, ycount = TRUE, min_male_ycount = 5) %>%
  dplyr::filter(status == "PROBLEM")
#> # A tibble: 3 × 6
#>   id    sex    snp_sex status        f y_count
#>   <chr> <fct>  <fct>   <chr>     <dbl>   <int>
#> 1 ind1  male   female  PROBLEM -0.0251       0
#> 2 ind2  female NA      PROBLEM  0.224        0
#> 3 ind20 female male    PROBLEM  1           10
```

## Analyses use autosomes only

Most analyses assume that loci are on diploid autosomes. On the X
chromosome, males are homozygous, which looks like inbreeding and adds
the sex of individuals to the structure of the data. As PLINK does
(e.g. for KING), the functions that pool loci into statistics of
individuals, pairs of individuals or populations, the PCA, and the
inputs of ADMIXTURE, sNMF and ADMIXTOOLS 2 use autosomal loci only, and
say how many loci they left out:

``` r

het <- indiv_het_obs(sex_gt)
#> Using autosomal loci only: 115 loci on X, Y, MT are left out (see ?gt_chromosome_types)
pca <- gt_pca(sex_gt, k = 2)
#> Using autosomal loci only: 115 loci on X, Y, MT are left out (see ?gt_chromosome_types)
```

The documentation of
[`gt_set_chromosome_types()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_chromosome_types.md)
lists these functions. To include a sex chromosome or mtDNA in them,
declare it as an autosome:

``` r

sex_gt %>%
  gt_set_chromosome_types(types = c(X = "autosome")) %>%
  gt_pca(k = 2)
#> Using autosomal loci only: 15 loci on Y, MT are left out (see ?gt_chromosome_types)
#>  === PCA of gen_tibble object ===
#> Method: [1] "partial"
#> 
#> Genotypes: complete 
#> 
#> Call ($call):gt_pca(x = ., k = 2)
#> 
#> Eigenvalues ($d):
#>  20.943 18.069 
#> 
#> Principal component scores ($u):
#>  matrix with 20 rows (individuals) and 2 columns (axes) 
#> 
#> Loadings (Principal axes) ($v):
#>  matrix with 150 rows (SNPs) and 2 columns (axes)
```

The functions that return one value per locus
(e.g. [`loci_alt_freq()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_alt_freq.md)
or
[`loci_hwe()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_hwe.md))
use all loci, but do not yet treat males as haploid: on X, they count
two alleles per male, and their HWE tests fail. Select the loci of the
types you need before using them.

## Writing sex chromosomes

[`gt_as_vcf()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_vcf.md)
writes the genotypes of males on X and Y (of females on Z and W), and of
everyone on mtDNA, as haploid genotypes. A heterozygous genotype of a
hemizygous individual (e.g. a genotyping error, or a male in the
pseudo-autosomal region if it was not declared) is written as missing,
as PLINK treats such genotypes as missing.

``` r

vcf_file <- gt_as_vcf(sex_gt, file = tempfile(fileext = ".vcf"))
#> 37 heterozygous genotypes of hemizygous individuals (e.g. males on X) were written as missing
vcf_lines <- readLines(vcf_file)
# the first genotypes of a locus on X (individuals 1 and 20 are recorded with
# the wrong sex)
x_line <- strsplit(vcf_lines[grep("^X\t", vcf_lines)[1]], "\t")[[1]]
x_line[10:29]
#>  [1] "1"   "0/1" "1/1" "1/1" "1/1" "1/1" "0/1" "0/1" "1/1" "0/1" "1"   "1"  
#> [13] "1"   "1"   "1"   "1"   "1"   "1"   "0"   "0/0"
```

[`gt_as_plink()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_as_plink.md)
writes the names of the chromosomes, or, with
`chromosomes_as_int = TRUE`, PLINK’s integer codes (23 to 26 for X, Y,
XY and MT in humans, or the codes of the number of autosomes declared
with `chr_set`).
