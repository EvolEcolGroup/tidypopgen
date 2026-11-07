# Write site frequency spectra for dadi or fastsimcoal2

Write site frequency spectra in the formats used by `dadi` (and
`moments`) and by `fastsimcoal2`, as written by `easySFS`:

- for an `sfs` object (as returned by
  [`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md),
  [`pairwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_pop_sfs.md)
  and
  [`nwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/nwise_pop_sfs.md)
  with `type = "array"`, or by the `bootstrap_*` functions, e.g.
  [`bootstrap_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/bootstrap_sfs.md)),
  write that spectrum to `file`;

- for a grouped
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md),
  compute all the spectra that `easySFS` writes (the 1D spectrum of each
  population, the joint spectrum of each pair of populations, and the
  joint spectrum of all populations, in a single pass over the
  genotypes), and write them in `out_dir`, with the file names of
  `easySFS`, optionally with bootstrap replicates. Missing genotypes are
  dealt with by projecting each population down to `n_proj` allele
  copies (see
  [`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md)).

## Usage

``` r
write_sfs(x, ...)

# S3 method for class 'sfs'
write_sfs(
  x,
  file,
  format = c("dadi", "fastsimcoal2", "fastsimcoal2_multidimensional"),
  fastest = c("last", "first"),
  total_length = NULL,
  round = FALSE,
  overwrite = FALSE,
  ...
)

# S3 method for class 'gen_tbl'
write_sfs(
  x,
  out_dir = NULL,
  format = c("dadi", "fastsimcoal2"),
  n_proj = NULL,
  prefix = NULL,
  folded = TRUE,
  ancestral = NULL,
  total_length = NULL,
  round = FALSE,
  fastest = c("last", "first"),
  n_bootstrap = 0,
  bootstrap_block_size = NULL,
  bootstrap_size_unit = c("snp", "cM", "bp"),
  max_cells = 1e+07,
  overwrite = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(x), 1),
  ...
)

# Default S3 method
write_sfs(x, ...)
```

## Arguments

- x:

  an `sfs` object, or a grouped
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  (as obtained by using
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)).

- ...:

  further arguments passed to the methods.

- file:

  for an `sfs` object, the name of the file to write.

- format:

  the format of the files: `"dadi"` or `"fastsimcoal2"`; for an `sfs`
  object, also `"fastsimcoal2_multidimensional"` (the multidimensional
  format of fastsimcoal2, whatever the number of populations). See
  details.

- fastest:

  the order of the cells in dadi files and in the multidimensional
  fastsimcoal2 format: "last" (the default) keeps the order of the
  populations, so that the number of allele copies in the last
  population changes fastest, as in `dadi` and `easySFS`; "first" makes
  the number of copies in the first population change fastest. See
  details.

- total_length:

  the total number of sites sequenced, including monomorphic ones, used
  to set the number of monomorphic sites (the first cell). If `NULL`
  (the default), the spectra are written as they are.

- round:

  boolean, whether to round the number of sites in each cell to the
  nearest integer (as `easySFS` with `--dtype int`). Defaults to
  `FALSE`.

- overwrite:

  boolean, whether to overwrite files that already exist. Defaults to
  `FALSE`.

- out_dir:

  for a `gen_tibble`, the directory in which to write the files, or
  `NULL` (the default) to use the location and name of the backing
  files, followed by `_dadi` or `_fastsimcoal2`.

- n_proj:

  the number of allele copies to project each population to: either a
  single value used for all populations, a vector with one value per
  population (in the order of the groups), or a vector named with the
  group labels. If `NULL` (the default), each population is projected to
  the number of allele copies that maximises its number of segregating
  sites (see
  [`pop_sfs_preview()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs_preview.md)).

- prefix:

  for a `gen_tibble` written in fastsimcoal2 format, the prefix of the
  names of the 2D and multidimensional files, or `NULL` (the default) to
  use the name of the backing files.

- folded:

  boolean, whether to fold the spectrum (the default, as the ancestral
  allele is not usually known). If `FALSE`, the spectrum is polarised
  with `ancestral`.

- ancestral:

  the ancestral allele of each locus, used to polarise unfolded spectra:
  either the name of a column of the loci table (see
  [`show_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/show_loci.md)),
  or a character vector with one allele per locus (in the order of the
  loci table). Loci whose ancestral allele is missing (`NA`) or is
  neither of their alleles are dropped. If `NULL` (the default), the
  reference allele is taken to be ancestral, as in `easySFS`. Ignored,
  with a warning, for folded spectra.

- n_bootstrap:

  for a `gen_tibble`, the number of bootstrap replicates to write (0,
  the default, for none).

- bootstrap_block_size, bootstrap_size_unit:

  for a `gen_tibble`, the size of the blocks of loci resampled for the
  bootstrap, and its unit ("snp", the default, "cM" or "bp"); see
  [`bootstrap_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/bootstrap_sfs.md).

- max_cells:

  for a `gen_tibble`, the maximum number of cells of the joint spectrum
  of all populations; an error is raised for larger spectra.

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html)

- block_size:

  maximum number of loci read at once.

## Value

the path of the file (for an `sfs` object) or the paths of the files
written (for a `gen_tibble`), invisibly.

## Details

**dadi format** (`format = "dadi"`): a line with the dimensions of the
spectrum (the number of allele copies of each population plus 1),
`folded` or `unfolded` and the population labels in quotes; a line with
the cells of the spectrum; and a line with the mask (1 for masked cells,
see
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md)).
In dadi files, the cells are listed with the number of copies in the
last population changing fastest. With `fastest = "first"`, the
populations are written in reverse order, so that the number of copies
in the first population changes fastest and the file is still a valid
dadi file.

**fastsimcoal2 format** (`format = "fastsimcoal2"`): for a 1D spectrum,
as a `_DAFpop0.obs` (unfolded) or `_MAFpop0.obs` (folded) file; for a 2D
spectrum, as a `_jointDAFpop1_0.obs` or `_jointMAFpop1_0.obs` file, a
matrix with the second population in the rows (labelled `d1_0`, `d1_1`,
...) and the first population in the columns (labelled `d0_0`, `d0_1`,
...); for a spectrum of three or more populations, as a multidimensional
spectrum (`_MSFS.obs` file). For an `sfs` object,
`format = "fastsimcoal2_multidimensional"` writes the multidimensional
format whatever the number of populations. In the multidimensional
format, the cells are listed in the order set by `fastest`: "last" (the
default) has the number of copies in the last population changing
fastest, as written by `easySFS`; "first" has the first population
changing fastest. Please check which order your version of
`fastsimcoal2` expects.

**Monomorphic sites.** As in `easySFS`, `total_length` gives the total
number of sites sequenced (including monomorphic sites that are not in
the `gen_tibble`): the first cell of the spectrum (sites with no copies
of the alternate or minor allele) is then set to `total_length` minus
the number of sites in all the other cells.

**The spectra of a `gen_tibble`** are written in `out_dir`, with the
file names of `easySFS` (`pop` stands for the population labels, `p` for
their projections, and `M` is replaced by `D` for unfolded spectra):

- dadi: `pop-p.sfs` for each population, `pop1-pop2.sfs` for each pair,
  and `pop1-pop2-...-popK.sfs` for all populations;

- fastsimcoal2: `pop_MAFpop0.obs` for each population,
  `prefix_jointMAFpopj_i.obs` for each pair (where `i` and `j` are the
  0-based indices of the two populations, so that the second population
  of the pair is in the rows), and `prefix_MSFS.obs` for all
  populations.

`easySFS` writes these files in the sub-directories `dadi` and
`fastsimcoal2` of its output directory: use, for example,
`out_dir = file.path(dir, "dadi")` with `format = "dadi"` and
`out_dir = file.path(dir, "fastsimcoal2")` with
`format = "fastsimcoal2"` to get the same layout. Populations are in the
order of the groups (see
[`dplyr::group_keys()`](https://dplyr.tidyverse.org/reference/group_data.html)):
to change their order, which matters for `dadi` and `fastsimcoal2` (and
is set with `--order` in `easySFS`), change the order of the levels of
the grouping variable (e.g. with
[`factor()`](https://rdrr.io/r/base/factor.html)).

With `n_bootstrap` replicates of a block bootstrap (see
[`bootstrap_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/bootstrap_sfs.md)),
the spectra of each replicate are also written, with the same file
names, in the sub-directories `bootstrap_1`, `bootstrap_2`, ... of
`out_dir`, numbered with as many digits as `n_bootstrap` (e.g.
`bootstrap_001` to `bootstrap_100`), which can be used to estimate
confidence intervals with `fastsimcoal2`. Replicates use the same
projections as the spectra of all loci, and `total_length` and `round`
are applied to each of them. To draw the same replicates in both
formats, set the same seed (with
[`set.seed()`](https://rdrr.io/r/base/Random.html)) before writing each
of them.

## See also

[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md),
[`pop_sfs_preview()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs_preview.md),
[`bootstrap_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/bootstrap_sfs.md)

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# write a single spectrum
sfs_2d <- example_gt %>%
  pairwise_pop_sfs(n_proj = 2, type = "array")
dadi_file <- tempfile(fileext = ".sfs")
write_sfs(sfs_2d$pop1.pop2, dadi_file)
readLines(dadi_file)
#> [1] "3 3 folded \"pop1\" \"pop2\""                                                                          
#> [2] "1.783333333333334 0.9333333333333336 0.3083333333333334 1.6 1.066666666666667 0 0.3083333333333334 0 0"
#> [3] "1 0 0 0 0 1 0 1 1"                                                                                     
fastsimcoal2_file <- tempfile(fileext = "_jointMAFpop1_0.obs")
write_sfs(sfs_2d$pop1.pop2, fastsimcoal2_file, format = "fastsimcoal2")
readLines(fastsimcoal2_file)
#> [1] "1 observation"                                 
#> [2] "\td0_0\td0_1\td0_2"                            
#> [3] "d1_0\t1.783333333333334 1.6 0.3083333333333334"
#> [4] "d1_1\t0.9333333333333336 1.066666666666667 0"  
#> [5] "d1_2\t0.3083333333333334 0 0"                  

# write all the spectra written by easySFS, in the same directories
out_dir <- file.path(tempdir(), "example_sfs")
dadi_files <- example_gt %>%
  write_sfs(out_dir = file.path(out_dir, "dadi"), n_proj = 2)
basename(dadi_files)
#> [1] "pop1-2.sfs"         "pop2-2.sfs"         "pop3-2.sfs"        
#> [4] "pop1-pop2.sfs"      "pop1-pop3.sfs"      "pop2-pop3.sfs"     
#> [7] "pop1-pop2-pop3.sfs"
fastsimcoal2_files <- example_gt %>%
  write_sfs(
    out_dir = file.path(out_dir, "fastsimcoal2"), format = "fastsimcoal2",
    n_proj = 2, prefix = "example"
  )
basename(fastsimcoal2_files)
#> [1] "pop1_MAFpop0.obs"           "pop2_MAFpop0.obs"          
#> [3] "pop3_MAFpop0.obs"           "example_jointMAFpop1_0.obs"
#> [5] "example_jointMAFpop2_0.obs" "example_jointMAFpop2_1.obs"
#> [7] "example_MSFS.obs"          

# with 3 bootstrap replicates, resampling blocks of 2 SNPs
set.seed(1)
fastsimcoal2_files <- example_gt %>%
  write_sfs(
    out_dir = file.path(out_dir, "fastsimcoal2"), format = "fastsimcoal2",
    n_proj = 2, prefix = "example", n_bootstrap = 3,
    bootstrap_block_size = 2, overwrite = TRUE
  )
list.files(file.path(out_dir, "fastsimcoal2"))
#>  [1] "bootstrap_1"                "bootstrap_2"               
#>  [3] "bootstrap_3"                "example_MSFS.obs"          
#>  [5] "example_jointMAFpop1_0.obs" "example_jointMAFpop2_0.obs"
#>  [7] "example_jointMAFpop2_1.obs" "pop1_MAFpop0.obs"          
#>  [9] "pop2_MAFpop0.obs"           "pop3_MAFpop0.obs"          
```
