# Save a gen_tibble

Save the tibble (and update the backing files). The `gen_tibble` object
is saved to a file with extension *.gt*, together with update its *.rds*
and *.bk* files. Note that multiple *.gt* files can be linked to the
same *.rds* and *.bk* files; generally, this occurs when we create
multiple subsets of the data. The *.gt* file then stores the information
on what subset of the full dataset we are interested in, whilst the
*.rds* and *.bk* file store the full dataset. To reload a `gen_tibble`,
you can pass the name of the *.gt* file with
[`gt_load()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_load.md).
The *.gt* file is compressed with gzip at its fastest level, as saving
the information on millions of loci with the default level of
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) is slow.

## Usage

``` r
gt_save(x, file_name = NULL, quiet = FALSE)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)

- file_name:

  the file name, including the full path. If it does not end with *.gt*,
  the extension will be added.

- quiet:

  boolean to suppress information about the files

## Value

the file name and path of the *.gt* file, together with the *.rds* and
*.bk* files

## Details

If the names of the loci are integers (`names_as_int = TRUE`), their
original names are in a sidecar file next to the *.bk* file
(*.loci.tsv*, see
[`gt_lookup_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md));
if the chromosomes are integers (`chromosomes_as_int = TRUE`), their
names are in another sidecar file (*.chromosomes.tsv*, see
[`gt_lookup_chromosomes()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)).
Keep any sidecar files with the other files (if you move the files, move
them too, so that
[`gt_load()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_load.md)
with `reattach_to` finds them).

## See also

[`gt_load()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_load.md)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# remove some individuals
example_gt <- example_gt %>% filter(id != "a")

# save filtered gen_tibble object
gt_save(example_gt, file_name = paste0(tempfile(), "_example_filtered"))
#> 
#> gen_tibble saved to /tmp/RtmpggvtNv/file1f3741bc6330_example_filtered.gt
#> using FBM RDS: /tmp/RtmpggvtNv/file1f374dfcbece.rds
#> with FBM backing file: /tmp/RtmpggvtNv/file1f374dfcbece.bk
#> make sure that you do NOT delete those files!
#> to reload the gen_tibble in another session, use:
#> gt_load('/tmp/RtmpggvtNv/file1f3741bc6330_example_filtered.gt')
#> [1] "/tmp/RtmpggvtNv/file1f3741bc6330_example_filtered.gt"
#> [2] "/tmp/RtmpggvtNv/file1f374dfcbece.rds"                
#> [3] "/tmp/RtmpggvtNv/file1f374dfcbece.bk"                 
```
