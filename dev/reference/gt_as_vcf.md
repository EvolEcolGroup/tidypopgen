# Convert a `gen_tibble` to a VCF

This function write a VCF from a `gen_tibble`.

## Usage

``` r
gt_as_vcf(x, file = NULL, chunk_size = NULL, overwrite = FALSE)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  with population coded as 'population'

- file:

  the .vcf file name with a path, or NULL (the default) to use the
  location of the backing files.

- chunk_size:

  the number of loci processed at a time. Automatically set if left to
  NULL

- overwrite:

  logical, should the file be overwritten if it already exists?

## Value

the path of the .vcf file

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Write a vcf file
example_gt %>% gt_as_vcf()
#> [1] "/tmp/Rtmpj4bjft/file21ae359394e8.vcf"
```
