# Convert a `gen_tibble` to a VCF

This function write a VCF from a `gen_tibble`. Phased genotypes (see
`phased` in
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md))
are written as phased (e.g. `0|1`), and unphased ones as unphased (e.g.
`0/1`). The genotypes are read and written in chunks of loci, so that
this works on datasets that do not fit in memory.

## Usage

``` r
gt_as_vcf(x, file = NULL, chunk_size = NULL, overwrite = FALSE)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  with population coded as 'population'

- file:

  the .vcf file name with a path, or NULL (the default) to use the
  location of the backing files.

- chunk_size:

  the number of loci processed at a time (also used for the number of
  chromosomes processed at a time when writing the contig lines of the
  header). If left to NULL, 100,000 loci are processed at a time.

- overwrite:

  logical, should the file be overwritten if it already exists?

## Value

the path of the .vcf file

## Details

On sex chromosomes and the mitochondrial genome (see
[`gt_chromosome_types()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_chromosome_types.md)),
hemizygous genotypes (stored as homozygotes, as in PLINK) are written as
haploid genotypes (e.g. `1` rather than `1/1`): those of individuals of
the heterogametic sex (males with X and Y, females with Z and W; from
the `sex` column) on X (or Z) outside the pseudo-autosomal regions and
on Y (or W), and those of all individuals on MT (and on Y or W if their
sex is unknown). Heterozygous genotypes of hemizygous individuals are
written as missing (`.`), as PLINK treats them as missing, with a
message giving their number, and the genotypes of the homogametic sex on
Y (or W) are written as missing.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Write a vcf file
example_gt %>% gt_as_vcf()
#> [1] "/tmp/Rtmp8cEXfI/file20367456f6dc.vcf"
```
