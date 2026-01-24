# Load a gen_tibble

Load a `gen_tibble` previously saved with
[`gt_save()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_save.md).
If the *.rds* and *.bk* files have not been moved, they should be found
automatically. If they were moved, use `reattach_to` to point to the
*.rds* file (the *.bk* file needs to be in the same directory as the
*.rds* file).

## Usage

``` r
gt_load(file = NULL, reattach_to = NULL)
```

## Arguments

- file:

  the file name, including the full path. If it does not end with *.gt*,
  the extension will be added.

- reattach_to:

  the file name, including the full path, of the *.rds* file if it was
  moved. It assumes that the *.bk* file is found in the same path. You
  should be able to leave this to NULL unless you have moved the files.

## Value

a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)

## See also

[`gt_save()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_save.md)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# remove some individuals
example_gt_filtered <- example_gt %>% filter(id != "a")

# save the filtered gen_tibble object
backing_files <- gt_save(example_gt_filtered,
  file_name = paste0(tempfile(), "_example_filtered")
)
#> 
#> gen_tibble saved to /tmp/RtmpKZXfxT/file218032d33f72_example_filtered.gt
#> using FBM RDS: /tmp/RtmpKZXfxT/file2180193543b3.rds
#> with FBM backing file: /tmp/RtmpKZXfxT/file2180193543b3.bk
#> make sure that you do NOT delete those files!
#> to reload the gen_tibble in another session, use:
#> gt_load('/tmp/RtmpKZXfxT/file218032d33f72_example_filtered.gt')

# backing_files[1] contains the name of the saved .gt file
backing_files[1]
#> [1] "/tmp/RtmpKZXfxT/file218032d33f72_example_filtered.gt"

# To load the saved gen_tibble object, use the path to the saved .gt file
reloaded_gt <- gt_load(backing_files[1])

# And we have loaded the gt without individual "a"
reloaded_gt
#> # A gen_tibble: 6 loci
#> # A tibble:     6 × 3
#>   id    population  genotypes
#>   <chr> <chr>      <vctr_SNP>
#> 1 b     pop1        [2,1,...]
#> 2 c     pop2        [2,.,...]
#> 3 d     pop2        [1,0,...]
#> 4 e     pop1        [1,2,...]
#> 5 f     pop3        [0,0,...]
#> 6 g     pop3        [0,1,...]
```
