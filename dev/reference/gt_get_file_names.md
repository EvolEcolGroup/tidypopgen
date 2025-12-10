# Get the names of files storing the genotypes of a `gen_tibble`

A function to return the names of the files used to store data in a
`gen_tibble`. Specifically, this returns the .rds file storing the big

## Usage

``` r
gt_get_file_names(x)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)

## Value

a character vector with the names and paths of the two files

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# To retrieve the names of and paths to the .bk and .rds files use:
gt_get_file_names(example_gt)
#> [1] "/tmp/Rtmp7A6Hoh/file2ab348bcb649.rds"
#> [2] "/tmp/Rtmp7A6Hoh/file2ab348bcb649.bk" 
```
