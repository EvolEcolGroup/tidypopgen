# Concatenate VCF files of the same samples into a `gen_tibble`

This function reads several *VCF* files with the same samples but
different loci (e.g. one file per chromosome, as distributed by many
projects) into a single
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
as `bcftools concat` would, followed by
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).
It needs no space on disk beyond the `gen_tibble`: each file is read
twice, once for its loci and samples, and once for its genotypes, which
are written straight into the backing file of the `gen_tibble`. To
combine files with different samples, use
[`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md).

## Usage

``` r
gt_concat_vcfs(
  files,
  backingfile = NULL,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  allow_duplicates = FALSE,
  names_as_int = TRUE,
  chromosomes_as_int = FALSE,
  n_cores = 1,
  quiet = FALSE
)
```

## Arguments

- files:

  a character vector with the paths to two or more *VCF* files (.vcf or
  .vcf.gz, which are read faster if compressed with `bgzip`), in the
  order in which their loci should be in the `gen_tibble`.

- backingfile:

  the path, including the file name without extension, for the backing
  files of the `gen_tibble` (a .bk file for the genotypes, a .rds file
  for the FBM, and a .gt file for the `gen_tibble`). If NULL, the path
  of the first file is used, with "\_concat" added to its name (and a
  version number if those files already exist, as in
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).

- valid_alleles:

  a vector of valid allele values; it defaults to 'A','T', 'C' and 'G'.

- missing_alleles:

  a vector of values in the *VCF* files that are considered as missing
  alleles ('.' is always missing).

- allow_duplicates:

  boolean: whether to allow duplicated loci (the same chromosome and
  position) or duplicated loci names in the `gen_tibble` (see
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).

- names_as_int:

  boolean: whether the loci names are integers (the index of each locus
  in the `gen_tibble`) rather than the *VCF* IDs, which saves memory
  with many loci (see
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).
  Note that it defaults to TRUE, unlike in
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  (where it defaults to FALSE): concatenating *VCF* files is mostly used
  for large datasets, in which the loci names take most of the memory of
  the `gen_tibble` and are rarely needed (loci are identified by their
  chromosome and position). Set it to FALSE to keep the *VCF* IDs. With
  integer names, the *VCF* IDs (with the chromosomes and positions) are
  written to a sidecar file next to the backing file, from which
  [`gt_lookup_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)
  and
  [`gt_restore_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)
  get them back.

- chromosomes_as_int:

  boolean: whether the chromosomes are integer codes rather than a
  factor, which saves memory with many chromosomes (e.g. RAD loci
  assembled de novo); the names are written to a sidecar file next to
  the backing file (see
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  and
  [`gt_lookup_chromosomes()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)).
  Default is FALSE.

- n_cores:

  the number of cores used to decompress bgzipped files and to parse the
  genotypes. Defaults to 1.

- quiet:

  boolean: whether to omit the report of the concatenation.

## Value

a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
with the data of all files (saved, as with
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).

## Details

As with
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
only biallelic SNPs are read from each file. The loci of the
`gen_tibble` are those of the first file, followed by those of the
second file, and so on (as in `bcftools concat`), so the files should be
given in the order of their chromosomes. The same chromosome can be
split across several files (e.g. in chunks), but duplicated loci (the
same chromosome and position) are an error, unless
`allow_duplicates = TRUE`.

All files must have the same samples, but they can be in a different
order in each file: the samples of the `gen_tibble` are in the order of
the first file, and the genotypes of the other files are reordered to
match. Sample names must be unique within each file.

As in
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
the ploidy of each individual is estimated from the first variants of
each file; the ploidy of an individual in the `gen_tibble` is its
largest ploidy across files (e.g. with a file of the X chromosome in
which males are haploid, males are still diploid). Files in which the
genotypes of an individual are all missing in the first variants are
ignored for its ploidy (unless they all are).

Unlike
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
by default (`names_as_int = TRUE`) the loci are named with integers
(their index in the `gen_tibble`), rather than with their *VCF* IDs, to
save memory with large datasets. With `names_as_int = FALSE`, loci take
their *VCF* IDs (a missing ID '.' gives a name made of the chromosome
and the position, as in
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).

## See also

[`gt_merge_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_merge_vcfs.md)
to merge *VCF* files with different samples.

## Examples

``` r
# two VCF files with the same individuals, one for each chromosome
example_gt <- load_example_gt("gen_tbl")
vcf_chr1 <- gt_as_vcf(select_loci(example_gt, 1:4),
  file = tempfile(fileext = ".vcf")
)
vcf_chr2 <- gt_as_vcf(select_loci(example_gt, 5:6),
  file = tempfile(fileext = ".vcf")
)
# concatenate them, keeping the VCF IDs as loci names (by default, unlike
# gen_tibble(), loci are named with integers to save memory)
concat_gt <- gt_concat_vcfs(c(vcf_chr1, vcf_chr2),
  names_as_int = FALSE,
  backingfile = tempfile("concat_")
)
#> reading the loci of /tmp/RtmpmlLt0o/file1f1d3b30f263.vcf
#> reading the loci of /tmp/RtmpmlLt0o/file1f1d4f0c57f0.vcf
#> concatenating 2 VCF files with 7 samples: 6 loci (biallelic SNPs)
#> reading the genotypes of /tmp/RtmpmlLt0o/file1f1d3b30f263.vcf
#> reading the genotypes of /tmp/RtmpmlLt0o/file1f1d4f0c57f0.vcf
#> 
#> gen_tibble saved to /tmp/RtmpmlLt0o/concat_1f1d76512bdd.gt
#> using FBM RDS: /tmp/RtmpmlLt0o/concat_1f1d76512bdd.rds
#> with FBM backing file: /tmp/RtmpmlLt0o/concat_1f1d76512bdd.bk
#> make sure that you do NOT delete those files!
#> to reload the gen_tibble in another session, use:
#> gt_load('/tmp/RtmpmlLt0o/concat_1f1d76512bdd.gt')
concat_gt
#> # A gen_tibble: 6 loci
#> # A tibble:     7 × 2
#>   id     genotypes
#>   <chr> <vctr_SNP>
#> 1 a      [1,1,...]
#> 2 b      [2,1,...]
#> 3 c      [2,.,...]
#> 4 d      [1,0,...]
#> 5 e      [1,2,...]
#> 6 f      [0,0,...]
#> 7 g      [0,1,...]
show_loci(concat_gt)
#> # A tibble: 6 × 7
#>   big_index name  chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr> <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         1 rs1   chr1              3            0 A          T         
#> 2         2 rs2   chr1              5            0 T          C         
#> 3         3 rs3   chr1             65            0 C          NA        
#> 4         4 rs4   chr1            343            0 G          C         
#> 5         5 rs5   chr2             23            0 C          G         
#> 6         6 rs6   chr2            456            0 T          A         
```
