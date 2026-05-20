# Combine two gen_tibbles

This function combined two
[gen_tibble](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)s.
By defaults, it subsets the loci and swaps ref and alt alleles to make
the two datasets compatible (this behaviour can be switched off with
`as_is`). The first object is used as a "reference" , and SNPs in the
other dataset will be flipped and/or alleles swapped as needed. SNPs
that have different alleles in the two datasets (i.e. triallelic) will
also be dropped. There are also options (NOT default) to attempt strand
flipping to match alleles (often needed in human datasets from different
SNP chips), and remove ambiguous alleles (C/G and A/T) where the correct
strand can not be guessed.

## Usage

``` r
# S3 method for class 'gen_tbl'
rbind(
  ...,
  as_is = FALSE,
  flip_strand = FALSE,
  use_position = FALSE,
  quiet = FALSE,
  backingfile = NULL
)
```

## Arguments

- ...:

  two
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  objects. Note that this function can not take more objects, `rbind`
  has to be done sequentially for large sets of objects.

- as_is:

  boolean determining whether the loci should be left as they are before
  merging. If FALSE (the defaults), `rbind` will attempt to subset and
  swap alleles as needed.

- flip_strand:

  boolean on whether strand flipping should be checked to match the two
  datasets. If this is set to TRUE, ambiguous SNPs (i.e. A/T and C/G)
  will also be removed. It defaults to FALSE

- use_position:

  boolean of whether a combination of chromosome and position should be
  used for matching SNPs. By default, `rbind` uses the locus name, so
  this is set to FALSE. When using 'use_position=TRUE', make sure
  chromosomes are coded in the same way in both `gen_tibbles` (a mix of
  e.g. 'chr1', '1' or 'chromosome1' can be the reasons if an
  unexpectedly large number variants are dropped when merging).

- quiet:

  boolean whether to omit reporting to screen

- backingfile:

  the path and prefix of the files used to store the merged data (it
  will be a .RDS to store the `bigSNP` object and a .bk file as its
  backing file for the FBM)

## Value

a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
with the merged data.

## Details

rbind differs from merging data with plink, which swaps the order of
allele1 and allele2 according to minor allele frequency when merging
datasets. rbind flips and/or swaps alleles according to the reference
dataset, not according to allele frequency.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Create a second gen_tibble to merge
test_indiv_meta <- data.frame(
  id = c("x", "y", "z"),
  population = c("pop1", "pop1", "pop2")
)
test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, 0, 0, 0),
  c(2, 2, 0, 0, 1, 1)
)
test_loci <- data.frame(
  name = paste0("rs", 1:6),
  chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
  position = as.integer(c(3, 5, 65, 343, 23, 456)),
  genetic_dist = as.double(rep(0, 6)),
  allele_ref = c("A", "T", "C", "G", "C", "T"),
  allele_alt = c("T", "C", NA, "C", "G", "A")
)

test_gt <- gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  valid_alleles = c("A", "T", "C", "G"),
  quiet = TRUE
)

# Merge the datasets using rbind
merged_gt <- rbind(ref = example_gt, target = test_gt, flip_strand = TRUE)
#> harmonising loci between two datasets
#> flip_strand =  TRUE  ; remove_ambiguous =  TRUE 
#> -----------------------------
#> dataset: reference 
#> number of SNPs: 6 reduced to 2 
#> ( 4 are ambiguous, of which 4  were removed)
#> -----------------------------
#> dataset: target 
#> number of SNPs: 6 reduced to 2 
#> ( 0 were flipped to match the reference set)
#> ( 4 are ambiguous, of which 4 were removed)
#> 
#> gen_tibble saved to /tmp/RtmpylcdcX/gt_merged_1fba3afaa4b3.gt
#> using FBM RDS: /tmp/RtmpylcdcX/gt_merged_1fba3afaa4b3.rds
#> with FBM backing file: /tmp/RtmpylcdcX/gt_merged_1fba3afaa4b3.bk
#> make sure that you do NOT delete those files!
#> to reload the gen_tibble in another session, use:
#> gt_load('/tmp/RtmpylcdcX/gt_merged_1fba3afaa4b3.gt')

merged_gt
#> # A gen_tibble: 2 loci
#> # A tibble:     10 × 3
#>    id    population  genotypes
#>    <chr> <chr>      <vctr_SNP>
#>  1 a     pop1        [1,0,...]
#>  2 b     pop1        [1,0,...]
#>  3 c     pop2        [.,0,...]
#>  4 d     pop2        [0,0,...]
#>  5 e     pop1        [2,0,...]
#>  6 f     pop3        [0,0,...]
#>  7 g     pop3        [1,1,...]
#>  8 x     pop1        [1,0,...]
#>  9 y     pop1        [1,0,...]
#> 10 z     pop2        [2,0,...]
```
