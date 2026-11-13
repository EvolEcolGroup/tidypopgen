# Merge several VCF files into a `gen_tibble`

This function reads several *VCF* files with different samples (e.g. one
per population, or per sequencing batch) into a single
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
as `bcftools merge` would, followed by
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).
It is much faster than reading each file with
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
and merging them with
[rbind()](https://evolecolgroup.github.io/tidypopgen/dev/reference/rbind.gen_tbl.md)
(one pair at a time), and needs no space on disk beyond the merged
`gen_tibble`: each file is read twice, once for its loci (to match them
across files) and once for its genotypes, which are written straight
into the merged backing file.

## Usage

``` r
gt_merge_vcfs(
  files,
  variants = c("intersection", "union"),
  missing_to_ref = FALSE,
  alt_mismatch = c("drop", "split"),
  ref_mismatch = c("error", "drop", "swap"),
  force_samples = FALSE,
  backingfile = NULL,
  valid_alleles = c("A", "T", "C", "G"),
  missing_alleles = c("0", "."),
  allow_duplicates = FALSE,
  names_as_int = TRUE,
  chromosomes_as_int = FALSE,
  n_cores = 1,
  phased = FALSE,
  quiet = FALSE
)
```

## Arguments

- files:

  a character vector with the paths to two or more *VCF* files (.vcf or
  .vcf.gz, which are read faster if compressed with `bgzip`).

- variants:

  which sites to keep: "intersection" (the default) keeps the sites
  found in all files; "union" keeps all sites (see details).

- missing_to_ref:

  boolean: with `variants = "union"`, whether the genotypes of the
  samples of a file at sites that are not in that file are homozygous
  for the reference allele (as with `bcftools merge -0`), rather than
  missing (the default).

- alt_mismatch:

  what to do with loci at the same position that have the same REF but a
  different ALT allele: "drop" (the default) drops the site, "split"
  keeps them as separate loci (see details).

- ref_mismatch:

  what to do with loci at the same position that have a different REF
  allele: "error" (the default) stops with an error, "drop" drops the
  site, "swap" merges loci whose alleles are swapped, swapping their
  genotypes (see details).

- force_samples:

  boolean: whether to allow the same sample name in different files,
  prefixing the duplicated names by the index of their file (see
  details). Defaults to FALSE.

- backingfile:

  the path, including the file name without extension, for the backing
  files of the merged `gen_tibble` (a .bk file for the genotypes, a .rds
  file for the FBM, and a .gt file for the `gen_tibble`). If NULL, the
  path of the first file is used, with "\_merged" added to its name (and
  a version number if those files already exist, as in
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).

- valid_alleles:

  a vector of valid allele values; it defaults to 'A','T', 'C' and 'G'.

- missing_alleles:

  a vector of values in the *VCF* files that are considered as missing
  alleles ('.' is always missing).

- allow_duplicates:

  boolean: whether to allow duplicated loci (the same chromosome and
  position) or duplicated loci names in the merged `gen_tibble` (see
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).

- names_as_int:

  boolean: whether the loci names are integers (the index of each locus
  in the merged `gen_tibble`) rather than the *VCF* IDs, which saves
  memory with many loci (see
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).
  Note that it defaults to TRUE, unlike in
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  (where it defaults to FALSE): merging *VCF* files is mostly used for
  large datasets, in which the loci names take most of the memory of the
  `gen_tibble` and are rarely needed (loci are identified by their
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

- phased:

  boolean: whether to store phased genotypes, giving a phased
  `gen_tibble` (see `phased` in
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).
  All the files must then be phased: any unphased genotype (with '/') or
  individual that is not diploid is an error. Default is FALSE (the
  phase is ignored, as in
  [`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).
  The alleles of phased genotypes are never swapped (as in
  [rbind()](https://evolecolgroup.github.io/tidypopgen/dev/reference/rbind.gen_tbl.md)),
  so `ref_mismatch = "swap"` is not available.

- quiet:

  boolean: whether to omit the report of the merge.

## Value

a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
with the merged data (saved, as with
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).

## Details

As with
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
only biallelic SNPs are read from each file. Loci are matched across
files by chromosome and position (as in `bcftools merge`), so
chromosomes need to be named in the same way in all files (a mix of e.g.
'chr1' and '1' would match no loci). A locus is merged with a site at
the same chromosome and position if their alleles agree: the same REF
and ALT alleles, or a missing allele (e.g. the ALT '.' of a monomorphic
locus), which takes the allele of the other files (a site with the same
alleles is preferred to one with a missing allele). Loci whose alleles
do not agree are handled as follows:

- loci with the same REF but a different ALT allele (e.g. A/C and A/G):
  `bcftools merge` would make a multiallelic site, which a `gen_tibble`
  can not store, so by default the site is dropped
  (`alt_mismatch = "drop"`); with `alt_mismatch = "split"` (as
  `bcftools merge -m none`), the loci are kept as separate loci at the
  same position (which needs `allow_duplicates = TRUE`), each with
  missing genotypes for the samples of the other files.

- loci with a different REF allele (e.g. A/C and C/A): the REF allele
  should be the base of the reference genome at that position, so
  `bcftools merge` stops with an error, as does this function by default
  (`ref_mismatch = "error"`, listing the first few mismatches). With
  `ref_mismatch = "drop"` such sites are dropped; with
  `ref_mismatch = "swap"`, loci whose alleles are swapped (A/C and C/A)
  are merged, swapping the genotypes of the later file (as
  [rbind()](https://evolecolgroup.github.io/tidypopgen/dev/reference/rbind.gen_tbl.md)
  does), and other mismatches are handled as with `alt_mismatch`.
  Swapping is needed e.g. for *VCF* files exported by PLINK 1.9, which
  by default set the REF allele to the major allele.

Strands are not flipped (all *VCF* files are on the forward strand).
Loci at the same position within a file are never merged with each other
(as in `bcftools merge`); if a file has several loci at the same
position (e.g. a multiallelic site split into biallelic records, as by
`bcftools norm -m-`), each is merged with the site of the other files
that has its alleles, if any, before the others are handled as above.

The merged `gen_tibble` has the samples of the first file, followed by
those of the second file, and so on. By default
(`variants = "intersection"`), it keeps the sites found in all files, as
[rbind()](https://evolecolgroup.github.io/tidypopgen/dev/reference/rbind.gen_tbl.md)
does; with `variants = "union"`, it keeps all sites (as
`bcftools merge`), and the samples of the files without a site get
missing genotypes there, or homozygous reference genotypes with
`missing_to_ref = TRUE` (as `bcftools merge -0`). Sites are ordered by
chromosome (in order of first appearance across the files), then by
position.

Unlike
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
by default (`names_as_int = TRUE`) the loci are named with integers
(their index in the merged `gen_tibble`), rather than with their *VCF*
IDs, to save memory with large datasets. With `names_as_int = FALSE`,
each locus takes its name (*VCF* ID) from the first file where it is
found (a missing ID '.' gives a name made of the chromosome and the
position, as in
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).

Sample names must be unique across files: if they are not, the function
stops, unless `force_samples = TRUE`, in which case the duplicated names
are prefixed by the index of their file, as
`bcftools merge --force-samples` does (e.g. sample 'S3' of the second
file becomes '2:S3').

As in
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
the ploidy of each individual is estimated from the first variants of
its file. The individuals of different files can have different
ploidies.

## See also

[`gt_concat_vcfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_concat_vcfs.md)
to combine *VCF* files with the same samples and different loci (e.g.
one file per chromosome);
[rbind()](https://evolecolgroup.github.io/tidypopgen/dev/reference/rbind.gen_tbl.md)
to merge two `gen_tibble` objects (e.g. read from files of other
formats).

## Examples

``` r
vcf_a <- system.file("extdata", "pop_a.vcf", package = "tidypopgen")
vcf_b <- system.file("extdata", "pop_b.vcf", package = "tidypopgen")
# the REF and ALT alleles of some loci are swapped between these two files
# (which a VCF from PLINK can have), so we swap them back to merge them;
# we also keep the VCF IDs as loci names (by default, unlike gen_tibble(),
# loci are named with integers to save memory)
merged_gt <- gt_merge_vcfs(c(vcf_a, vcf_b),
  ref_mismatch = "swap",
  names_as_int = FALSE,
  backingfile = tempfile("merged_")
)
#> reading the loci of /home/runner/work/_temp/Library/tidypopgen/extdata/pop_a.vcf
#> reading the loci of /home/runner/work/_temp/Library/tidypopgen/extdata/pop_b.vcf
#> merging 2 VCF files with 8 samples:
#>   20 sites (biallelic SNPs), of which 13 are in all files
#>   dropped: 5 with alleles that would make a multiallelic site, 0 with different REF alleles
#>   loci with swapped alleles: 2
#>   kept: 8 loci (intersection)
#> reading the genotypes of /home/runner/work/_temp/Library/tidypopgen/extdata/pop_a.vcf
#> reading the genotypes of /home/runner/work/_temp/Library/tidypopgen/extdata/pop_b.vcf
#> 
#> gen_tibble saved to /tmp/RtmpJMtr74/merged_1fc26dd7523f.gt
#> using FBM RDS: /tmp/RtmpJMtr74/merged_1fc26dd7523f.rds
#> with FBM backing file: /tmp/RtmpJMtr74/merged_1fc26dd7523f.bk
#> make sure that you do NOT delete those files!
#> to reload the gen_tibble in another session, use:
#> gt_load('/tmp/RtmpJMtr74/merged_1fc26dd7523f.gt')
merged_gt
#> # A gen_tibble: 8 loci
#> # A tibble:     8 × 2
#>   id                 genotypes
#>   <chr>             <vctr_SNP>
#> 1 pop_a_GRC14300079  [0,0,...]
#> 2 pop_a_GRC14300142  [0,0,...]
#> 3 pop_a_GRC14300159  [1,1,...]
#> 4 pop_a_GRC14300286  [0,0,...]
#> 5 pop_a_GRC14300305  [0,1,...]
#> 6 pop_b_SA073        [0,0,...]
#> 7 pop_b_SA1021       [1,1,...]
#> 8 pop_b_SA1008       [0,1,...]
show_loci(merged_gt)
#> # A tibble: 8 × 7
#>   big_index name       chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr>      <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         1 rs3094315  1            752566            0 A          G         
#> 2         2 rs3131972  1            752721            0 G          A         
#> 3         3 rs12124819 1            776546            0 A          G         
#> 4         4 rs11240777 1            798959            0 G          A         
#> 5         5 rs1110052  1            873558            0 G          T         
#> 6         6 rs6657048  1            957640            0 C          T         
#> 7         7 rs2488991  1            994391            0 T          G         
#> 8         8 rs1240719  1           1346905            0 T          A         
```
