# Generate a report of what would happen to each SNP in a merge

This function provides an overview of the fate of each SNP in two
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
objects in the case of a merge. Only SNPs found in both objects will be
kept. One object is used as a `reference`, and SNPs in the other dataset
will be flipped and/or alleles swapped as needed. SNPs that have
different alleles in the two datasets will also be dropped.

## Usage

``` r
rbind_dry_run(
  ref,
  target,
  use_position = FALSE,
  flip_strand = FALSE,
  quiet = FALSE
)
```

## Arguments

- ref:

  either a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object, or the path to the PLINK bim file; the alleles in this objects
  will be used as template to flip the ones in `target` and/or swap
  their order as necessary.

- target:

  either a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object, or the path to the PLINK bim file

- use_position:

  boolean of whether a combination of chromosome and position should be
  used for matching SNPs. By default, `rbind` uses the locus name, so
  this is set to FALSE. When using 'use_position=TRUE', make sure
  chromosomes are coded in the same way in both `gen_tibbles` (a mix of
  e.g. 'chr1', '1' or 'chromosome1' can be the reasons if an
  unexpectedly large number variants are dropped when merging).

- flip_strand:

  boolean on whether strand flipping should be checked to match the two
  datasets. Ambiguous SNPs (i.e. A/T and C/G) will also be removed. It
  defaults to FALSE

- quiet:

  boolean whether to omit reporting to screen

## Value

a list with two `data.frames`, named `target` and `ref`. Each data.frame
has [`nrow()`](https://rdrr.io/r/base/nrow.html) equal to the number of
loci in the respective dataset, a column `id` with the locus name, and
boolean columns `to_keep` (the valid loci that will be kept in the
merge), `alleles_mismatched` (loci found in both datasets but with
mismatched alleles, leading to those loci being dropped), `to_flip`
(loci that need to be flipped to align the two datasets, only found in
`target` data.frame) and `to_swap` (loci for which the order of alleles
needs to be swapped to align the two datasets, `target` data.frame)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Create a second gen_tibble to merge
test_indiv_meta <- data.frame(
  id = c("x", "y", "z"),
  population = c("pop1", "pop1", "pop2")
)
test_genotypes <- rbind(
  c(1, 1, 2, 1, 1),
  c(2, 1, 2, 0, 0),
  c(2, 2, 2, 0, 1)
)
test_loci <- data.frame(
  name = paste0("rs", 1:5),
  chromosome = paste0("chr", c(1, 1, 1, 1, 2)),
  position = as.integer(c(3, 5, 65, 343, 23)),
  genetic_dist = as.double(rep(0, 5)),
  allele_ref = c("A", "T", "C", "G", "C"),
  allele_alt = c("T", "C", NA, "C", "G")
)

test_gt <- gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  valid_alleles = c("A", "T", "C", "G"),
  quiet = TRUE
)

# Create an rbind report using rbind_dry_run
rbind_dry_run(example_gt, test_gt, flip_strand = TRUE)
#> harmonising loci between two datasets
#> flip_strand =  TRUE  ; remove_ambiguous =  TRUE 
#> -----------------------------
#> dataset: reference 
#> number of SNPs: 6 reduced to 2 
#> ( 4 are ambiguous, of which 4  were removed)
#> -----------------------------
#> dataset: target 
#> number of SNPs: 5 reduced to 2 
#> ( 0 were flipped to match the reference set)
#> ( 3 are ambiguous, of which 3 were removed)
#> $target
#>   id new_id name to_flip to_swap missing_allele ambiguous
#> 1  1     NA  rs1   FALSE   FALSE             NA      TRUE
#> 2  2      1  rs2   FALSE   FALSE             NA     FALSE
#> 3  3      2  rs3   FALSE   FALSE             NA     FALSE
#> 4  4     NA  rs4   FALSE   FALSE             NA      TRUE
#> 5  5     NA  rs5   FALSE   FALSE             NA      TRUE
#> 
#> $ref
#>   id new_id name missing_allele ambiguous
#> 1  1     NA  rs1             NA      TRUE
#> 2  2      1  rs2             NA     FALSE
#> 3  3      2  rs3             NA     FALSE
#> 4  4     NA  rs4             NA      TRUE
#> 5  5     NA  rs5             NA      TRUE
#> 6  6     NA  rs6             NA      TRUE
#> 
#> attr(,"class")
#> [1] "rbind_report" "list"        
#> attr(,"flip_strand")
#> [1] TRUE
#> attr(,"remove_ambiguous")
#> [1] TRUE
```
