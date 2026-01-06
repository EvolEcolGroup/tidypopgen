# Convert a `gentibble` to a .geno file for sNMF from the LEA package

This function writes a .geno file from a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).
Unless a file path is given, a file with suffix .geno is written in the
same location as the .rds and .bk files that underpin the
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).

## Usage

``` r
gt_as_geno_lea(x, file = NULL)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)

- file:

  the .geno filename with a path, or NULL (the default) to use the
  location of the backing files.

## Value

the path of the .geno file

## Details

NOTE that we currently read all the data into memory to write the file,
so this function is not suitable for very large datasets.

## See also

[`LEA::geno()`](https://rdrr.io/pkg/LEA/man/geno.html)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Write a geno file
gt_as_geno_lea(example_gt, file = paste0(tempfile(), "_example.geno"))
#> [1] "/tmp/RtmpQ5xhg7/file21934dcc4548_example.geno"
```
