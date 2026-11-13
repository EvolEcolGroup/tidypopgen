# Get the names of files storing the genotypes of a `gen_tibble`

A function to return the names of the files used to store data in a
`gen_tibble`. Specifically, this returns the *.rds* file storing the big
matrix (FBM) and its *.bk* backing file, always as the first two
elements and in this order. If the `gen_tibble` has integer names of
loci (`names_as_int = TRUE`) or integer chromosomes
(`chromosomes_as_int = TRUE`), the sidecar files holding their original
names (*.loci.tsv*, then *.chromosomes.tsv*, see
[`gt_lookup_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)
and
[`gt_lookup_chromosomes()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md))
follow. A sidecar is listed only if the `gen_tibble` uses it (its loci
table holds integer names, or integer chromosomes) and the file exists,
so a `gen_tibble` without integer names or chromosomes has exactly two
files.

## Usage

``` r
gt_get_file_names(x)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)

## Value

a character vector with the names and paths of the *.rds* and *.bk*
files, followed by those of any sidecar files

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# To retrieve the names of and paths to the .bk and .rds files use:
gt_get_file_names(example_gt)
#> [1] "/tmp/RtmpJMtr74/file1fc26fa3f4e1.rds"
#> [2] "/tmp/RtmpJMtr74/file1fc26fa3f4e1.bk" 
```
