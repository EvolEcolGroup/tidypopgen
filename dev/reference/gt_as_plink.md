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
#> [1] "/tmp/RtmpkB239F/file1f0555e45887_plink.bed"

# Write a ped file
example_gt %>% gt_as_plink(type = "ped", file = paste0(tempfile(), "_plink"))
#> [1] "/tmp/RtmpkB239F/file1f056fa7e0c5_plink.ped"

# Write a raw file
example_gt %>% gt_as_plink(type = "raw", file = paste0(tempfile(), "_plink"))
#> [1] "/tmp/RtmpkB239F/file1f0526a313e4_plink.raw"
```
