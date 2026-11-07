# Convert a `gentibble` to a .geno file for sNMF from the LEA package

This function writes a .geno file from a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).
Unless a file path is given, a file with suffix .geno is written in the
same location as the .rds and .bk files that underpin the
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

## Usage

``` r
gt_as_geno_lea(x, file = NULL)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)

- file:

  the .geno filename with a path, or NULL (the default) to use the
  location of the backing files.

## Value

the path of the .geno file

## Details

The genotypes are read and written directly from the file backed matrix
(one locus at a time), so this function also works on datasets that do
not fit in memory.

## Note

the `.geno` format only supports haploid or diploid data, so this
function only works on diploid (or pseudohaploid) `gen_tibble` objects.
Pseudohaploid data is simply treated as regular diploid genotype counts
(dosage 0/2), the same convention used by other tools (e.g. PLINK).

## See also

[`LEA::geno()`](https://rdrr.io/pkg/LEA/man/geno.html)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Write a geno file
gt_as_geno_lea(example_gt, file = paste0(tempfile(), "_example.geno"))
#> [1] "/tmp/RtmpRxbawj/file1ecf5934ac3d_example.geno"
```
