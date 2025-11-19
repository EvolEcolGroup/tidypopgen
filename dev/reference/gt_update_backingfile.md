# Update the backing matrix

This functions forces a re-write of the file backing matrix to match the
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).
Individuals and loci are subsetted and reordered according to the
current state of the `gen_tibble`. Tests for this function are in
test_gt_order_loci.R

## Usage

``` r
gt_update_backingfile(
  .x,
  backingfile = NULL,
  chunk_size = NULL,
  rm_unsorted_dist = TRUE,
  quiet = FALSE
)
```

## Arguments

- .x:

  a `gen_tibble` object

- backingfile:

  the path, including the file name without extension, for backing files
  used to store the data (they will be given a .bk and .RDS
  automatically). If left to NULL (the default), the file name will be
  based on the name f the current backing file.

- chunk_size:

  the number of loci to process at once

- rm_unsorted_dist:

  boolean to set `genetic_dist` to zero (i.e. remove it) if it is
  unsorted within the chromosomes.

- quiet:

  boolean to suppress information about the files

## Value

a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
with a backing file (i.e. a new File Backed Matrix)

## Details

This function does not check whether the positions of your genetic loci
are sorted. To check this, and update the file backing matrix, use
[`gt_order_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_order_loci.md).

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% gt_update_backingfile()
#> 
#> gen_backing files updated, now
#> using FBM RDS: /tmp/Rtmpyikrxm/file1f10132bcfeb_v2.rds
#> with FBM backing file: /tmp/Rtmpyikrxm/file1f10132bcfeb_v2.bk
#> make sure that you do NOT delete those files!
#> # A gen_tibble: 6 loci
#> # A tibble:     7 × 3
#>   id    population  genotypes
#>   <chr> <chr>      <vctr_SNP>
#> 1 a     pop1        [1,1,...]
#> 2 b     pop1        [2,1,...]
#> 3 c     pop2        [2,.,...]
#> 4 d     pop2        [1,0,...]
#> 5 e     pop1        [1,2,...]
#> 6 f     pop3        [0,0,...]
#> 7 g     pop3        [0,1,...]
```
