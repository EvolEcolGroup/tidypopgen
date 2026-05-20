# Update the backing matrix

This functions forces a re-write of the file backing matrix to match the
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).
Individuals and loci are subsetted and reordered according to the
current state of the `gen_tibble`. A `.gt`, `.bk` and `.rds` file will
be created.

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
  used to store the data (they will be given a .bk and .rds
  automatically). If left to NULL (the default), the file name will be
  based on the name of the current backing file.

- chunk_size:

  the number of loci to process at once

- rm_unsorted_dist:

  boolean to set `genetic_dist` to zero (i.e. remove it) if it is
  unsorted within the chromosomes.

- quiet:

  boolean to suppress information about the files

## Value

a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
with updated `.gt`, `.bk`, and `.rds` files (i.e. a new File Backed
Matrix)

## Details

This function does not check whether the positions of your genetic loci
are sorted. To check this, and update the file backing matrix, use
[`gt_order_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_order_loci.md).
Tests for this function are in test_gt_order_loci.R

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Here, new backingfiles are created, but without using `<-` the object
# loaded in the R session is not updated
gt_update_backingfile(example_gt)
#> 
#> gen_backing files updated, now
#> using FBM RDS: /tmp/RtmpylcdcX/file1fba52cb924a_v2.rds
#> with FBM backing file: /tmp/RtmpylcdcX/file1fba52cb924a_v2.bk
#> make sure that you do NOT delete those files!
#> to reload the gen_tibble in another session, use:
#> gt_load('/tmp/RtmpylcdcX/file1fba52cb924a_v2.gt')
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

# Make sure to use `<-` to update the file names associated with the
# gen_tibble object loaded in the R session
example_gt <- example_gt %>% gt_update_backingfile()
#> 
#> gen_backing files updated, now
#> using FBM RDS: /tmp/RtmpylcdcX/file1fba52cb924a_v3.rds
#> with FBM backing file: /tmp/RtmpylcdcX/file1fba52cb924a_v3.bk
#> make sure that you do NOT delete those files!
#> to reload the gen_tibble in another session, use:
#> gt_load('/tmp/RtmpylcdcX/file1fba52cb924a_v3.gt')
```
