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
resulting plink files. NOTE that writing to bed has been optimised for
speed, but writing to ped or raw is slower, especially for large
datasets.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Write a bed file
example_gt %>% gt_as_plink(type = "bed", file = paste0(tempfile(), "_plink"))
#> [1] "/tmp/RtmpNWzAb2/file21562d96dfdf_plink.bed"

# Write a ped file
example_gt %>% gt_as_plink(type = "ped", file = paste0(tempfile(), "_plink"))
#> [1] "/tmp/RtmpNWzAb2/file2156483f0626_plink.ped"

# Write a raw file
example_gt %>% gt_as_plink(type = "raw", file = paste0(tempfile(), "_plink"))
#> [1] "/tmp/RtmpNWzAb2/file215621e481c0_plink.raw"
```
