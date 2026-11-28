# Export a `gen_tibble` object to PLINK bed format

This function exports all the information of a `gen_tibble` object into
a PLINK bed, ped or raw file (and associated files, i.e. .bim and .fam
for .bed; .fam for .ped).

## Usage

``` r
gt_as_plink(
  x,
  file = NULL,
  type = c("bed", "ped", "raw"),
  overwrite = TRUE,
  chromosomes_as_int = FALSE
)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object

- file:

  a character string giving the path to output file. If left to NULL,
  the output file will have the same path and prefix of the backingfile.

- type:

  one of "bed", "ped" or "raw"

- overwrite:

  boolean whether to overwrite the file.

- chromosomes_as_int:

  boolean whether to use the integer representation of the chromosomes
  (for "bed" files): numbers are kept, and sex chromosomes and mtDNA
  (see
  [`gt_chromosome_types()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_chromosome_types.md))
  get PLINK's codes: with n autosomes, X (or Z) is n+1, Y (or W) n+2,
  the pseudo-autosomal regions (XY) n+3 and MT n+4, where n is the
  number of autosomes declared with `chr_set` in
  [`gt_set_chromosome_types()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_chromosome_types.md),
  or 22 (PLINK's default, for humans), or the largest number of an
  autosome if it is larger; other names are numbered after n+4.

## Value

the path of the saved file

## Details

If the gen_tibble has been read in from vcf format, family.ID in the
resulting plink files will be the same as sample.ID. If the gen_tibble
has a grouping variable, this will be used as the family.ID in the
resulting plink files (with several grouping variables, their values are
pasted with ".", e.g. `pop1.F`); otherwise, the population is used as
the family.ID. Parents, sex and phenotype are taken from the columns
`paternal_ID`, `maternal_ID`, `sex` and `phenotype`, if present, as
created by
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
when reading PLINK files (bed, ped or bigSNP rds): sex is written as 1
for "male", 2 for "female" and 0 (unknown) otherwise; phenotype as 1 for
"control" and 2 for "case", and other values (e.g. a quantitative
phenotype) as they are; missing parents are written as 0, and missing
phenotypes as -9. Genotypes are read and written directly from the file
backed matrix, so all formats can be written for datasets that do not
fit in memory (but note that ped and raw files are much larger, and thus
slower to write, than bed files).

## Note

PLINK's bed, ped and raw formats only support diploid data, so this
function only works on diploid (or pseudohaploid) `gen_tibble` objects.
Pseudohaploid data is simply treated as regular diploid genotype counts
(dosage 0/2), the same convention used by PLINK itself.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Write a bed file
example_gt %>% gt_as_plink(type = "bed", file = paste0(tempfile(), "_plink"))
#> [1] "/tmp/Rtmpo1kR2l/file1f4d3827ecbe_plink.bed"

# Write a ped file
example_gt %>% gt_as_plink(type = "ped", file = paste0(tempfile(), "_plink"))
#> [1] "/tmp/Rtmpo1kR2l/file1f4d63e1b737_plink.ped"

# Write a raw file
example_gt %>% gt_as_plink(type = "raw", file = paste0(tempfile(), "_plink"))
#> [1] "/tmp/Rtmpo1kR2l/file1f4d7b68003_plink.raw"
```
