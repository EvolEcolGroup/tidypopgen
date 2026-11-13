# Look up the original names of loci and chromosomes stored as integers

With `names_as_int = TRUE` and `chromosomes_as_int = TRUE`,
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
stores the names of the loci, and their chromosomes, as integers, to
save memory and time with large datasets (e.g. RAD loci assembled de
novo, with millions of loci, each its own "chromosome"). The original
names are written to sidecar files next to the backing file (see
Details), from which these functions get them back:

- `gt_lookup_loci()` gives back the original names and chromosomes of
  some loci (e.g. the hits of a selection scan), without changing the
  `gen_tibble`;

- `gt_lookup_chromosomes()` gives back the names of chromosome codes
  (e.g. the `chromosome` column of the results of window statistics or
  runs of homozygosity);

- `gt_restore_loci()` puts the original names back into the loci table
  of a `gen_tibble` (e.g. once it has been filtered down to fewer loci).

## Usage

``` r
gt_lookup_loci(.x, loci, check_positions = TRUE)

gt_lookup_chromosomes(.x, chromosomes)

gt_restore_loci(.x, names = TRUE, chromosomes = TRUE, check_positions = TRUE)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  (or, except for `gt_restore_loci()`, its `genotypes` column)

- loci:

  a subset of the loci table of `.x`: rows of
  [`show_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/show_loci.md),
  in any order, possibly with extra columns (e.g. the statistics of a
  scan); or row numbers of `show_loci(.x)`, as given by
  [`which()`](https://rdrr.io/r/base/which.html).

- check_positions:

  boolean: whether to check that the positions of the loci in the locus
  sidecar are the same as in `loci` (or in the loci table of `.x`, for
  `gt_restore_loci()`). For `gt_lookup_loci()`, they are only checked if
  `loci` has a `position` column. Default is TRUE.

- chromosomes:

  for `gt_lookup_chromosomes()`, a vector of chromosome codes (e.g. the
  `chromosome` column of the results of
  [`windows_pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_pairwise_pop_fst.md)
  or
  [`gt_roh_window()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_roh_window.md));
  for `gt_restore_loci()`, boolean: whether to restore the names of the
  chromosomes (if they are integer codes), default TRUE.

- names:

  boolean: whether to restore the names of the loci (if they are
  integers). Default is TRUE.

## Value

`gt_lookup_loci()` returns `loci` (as a table), with integer names in
its `name` column, and integer codes in its `chromosome` column,
replaced by their original names (as character vectors); the other
columns, and the rows and their order, are unchanged. Columns that
already hold names are left as they are, so the same code works on data
read with or without the integer options.

`gt_lookup_chromosomes()` returns the names of the chromosomes (a
character vector, with NA for NA). Chromosomes that are already names (a
factor or a character vector) are returned as a character vector.

`gt_restore_loci()` returns `.x`, with the original names in its loci
table (as if it had been read without the integer options: names as a
character vector, and chromosomes as a factor). The levels of the factor
are the chromosomes of the loci in `.x` (unused levels are dropped, as
by [`droplevels()`](https://rdrr.io/r/base/droplevels.html)), so that
restoring a few loci of data with millions of chromosomes is quick. It
does not save it: use
[`gt_save()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_save.md).

## Details

The sidecars of the backing file `foo.bk` are:

- `foo.loci.tsv` (written with `names_as_int = TRUE`), with the columns
  `int_name`, `name`, `chromosome` and `position`, and one line per
  locus, in increasing order of integer name. The names are those that
  `names_as_int = FALSE` would give (with a missing ID replaced by the
  chromosome and position), and the chromosomes are their names (not
  their codes);

- `foo.chromosomes.tsv` (written with `chromosomes_as_int = TRUE`), with
  the columns `int_chromosome` and `chromosome`, and one line per
  chromosome, in increasing order of code.

Both are tab separated text files with a header, with no quotes, and "."
for missing values, so they can be read with any tool (e.g.
`data.table::fread(path, quote = "")` or
`readr::read_tsv(path, quote = "")`), and joined to
[`show_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/show_loci.md)
by `name = int_name` or `chromosome = int_chromosome`.

The sidecars are written by
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
when it reads a VCF file with `parser = "cpp"`, and by
[`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md)
and
[`gt_concat_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_concat_vcfs.md).
Functions that write a new backing file
([`gt_update_backingfile()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_update_backingfile.md),
[`gt_order_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_order_loci.md),
[rbind()](https://evolecolgroup.github.io/tidypopgen/dev/reference/rbind.gen_tbl.md))
write new sidecars for it, with the loci and chromosomes that they keep
(integer names and codes keep their values, so they still find their
line). Keep the sidecars with the other files of the `gen_tibble` (see
[`gt_save()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_save.md)):
without them, the original names are lost.

Only the file name ties a sidecar to its data, so the position of each
locus in the locus sidecar is compared with its position in the loci
table: a different position means that the sidecar belongs to other data
(e.g. files copied by hand from different runs), and is an error. If you
changed the positions of the loci yourself (e.g. after a liftover to a
new assembly), the names in the sidecar are still right: use
`check_positions = FALSE`. The chromosome sidecar has no positions to
check.

Lookups are fast: the first lookup in a session reads through a sidecar
once to index it (about a second per GB), and later lookups of a few
loci only read the parts of the file that hold them.

## Examples

``` r
vcf_path <- system.file(
  "extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
  package = "tidypopgen"
)
anole_gt <- gen_tibble(vcf_path,
  backingfile = tempfile("anolis_"),
  names_as_int = TRUE, chromosomes_as_int = TRUE, quiet = TRUE
)
head(show_loci(anole_gt))
#> # A tibble: 6 × 7
#>   big_index  name chromosome position genetic_dist allele_ref allele_alt
#>       <int> <int>      <int>    <int>        <dbl> <chr>      <chr>     
#> 1         1     1          1        2            0 G          A         
#> 2         2     2          1       13            0 A          G         
#> 3         3     3          1       16            0 T          C         
#> 4         4     4          2       57            0 G          A         
#> 5         5     5          3       34            0 T          C         
#> 6         6     6          3       37            0 A          C         

# the original names of the loci with the most missing genotypes
missing <- loci_missingness(anole_gt)
gt_lookup_loci(anole_gt, which(missing > 0.6))
#> # A tibble: 81 × 7
#>    big_index name         chromosome position genetic_dist allele_ref allele_alt
#>        <int> <chr>        <chr>         <int>        <dbl> <chr>      <chr>     
#>  1        25 locus_2327_… locus_2327       59            0 A          C         
#>  2        71 locus_6428_… locus_6428       40            0 C          T         
#>  3        73 locus_6428_… locus_6428       42            0 T          G         
#>  4       236 locus_19280… locus_192…        4            0 G          T         
#>  5       254 locus_19654… locus_196…       88            0 C          A         
#>  6       255 locus_19654… locus_196…       89            0 C          A         
#>  7       256 locus_19654… locus_196…       90            0 A          T         
#>  8       349 locus_21964… locus_219…       78            0 A          G         
#>  9       364 locus_22697… locus_226…      128            0 C          T         
#> 10       451 locus_26111… locus_261…       21            0 G          T         
#> # ℹ 71 more rows

# the same, from a subset of the loci table, with the missingness
hits <- show_loci(anole_gt) %>%
  dplyr::mutate(missing = missing) %>%
  dplyr::filter(missing > 0.6)
gt_lookup_loci(anole_gt, hits)
#> # A tibble: 81 × 8
#>    big_index name         chromosome position genetic_dist allele_ref allele_alt
#>        <int> <chr>        <chr>         <int>        <dbl> <chr>      <chr>     
#>  1        25 locus_2327_… locus_2327       59            0 A          C         
#>  2        71 locus_6428_… locus_6428       40            0 C          T         
#>  3        73 locus_6428_… locus_6428       42            0 T          G         
#>  4       236 locus_19280… locus_192…        4            0 G          T         
#>  5       254 locus_19654… locus_196…       88            0 C          A         
#>  6       255 locus_19654… locus_196…       89            0 C          A         
#>  7       256 locus_19654… locus_196…       90            0 A          T         
#>  8       349 locus_21964… locus_219…       78            0 A          G         
#>  9       364 locus_22697… locus_226…      128            0 C          T         
#> 10       451 locus_26111… locus_261…       21            0 G          T         
#> # ℹ 71 more rows
#> # ℹ 1 more variable: missing <dbl>

# the names of some chromosome codes
gt_lookup_chromosomes(anole_gt, c(1, 5, 5, NA))
#> [1] "locus_137" "locus_575" "locus_575" NA         

# put the original names back into the gen_tibble
anole_gt <- gt_restore_loci(anole_gt)
head(show_loci(anole_gt))
#> # A tibble: 6 × 7
#>   big_index name         chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr>        <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         1 locus_137_2  locus_137         2            0 G          A         
#> 2         2 locus_137_13 locus_137        13            0 A          G         
#> 3         3 locus_137_16 locus_137        16            0 T          C         
#> 4         4 locus_203_57 locus_203        57            0 G          A         
#> 5         5 locus_346_34 locus_346        34            0 T          C         
#> 6         6 locus_346_37 locus_346        37            0 A          C         
```
