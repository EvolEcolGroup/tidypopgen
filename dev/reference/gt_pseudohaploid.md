# Set the ploidy of a `gen_tibble` to include pseudohaploids

If a `gen_tibble` includes pseudohaploid data, its ploidy is set to -2
to indicate that some individuals are coded as pseudohaploids. The
ploidy of the individuals is updated, with pseudohaploids set to 1 and
diploids set to 2. However, the dosages are not changed, meaning that
pseudohaploids are still coded as 0 or 2. If the `gen_tibble` is already
set to pseudohaploid, running gt_pseudohaploid will update the ploidy
values again, if pseudohaploid individuals have been removed then ploidy
is reset to 2.

## Usage

``` r
gt_pseudohaploid(x, test_n_loci = 10000)
```

## Arguments

- x:

  a `gen_tibble` object

- test_n_loci:

  the number of loci to test to determine if an individual is
  pseudohaploid. If there are no heterozygotes in the first
  `test_n_loci` loci, the individual is considered a pseudohaploid. If
  `NULL`, all loci are tested.

## Value

a `gen_tibble` object with the ploidy set to -2 and the individual
ploidy values updated to 1 or 2.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Detect pseudohaploids and set ploidy for the whole gen_tibble
example_gt <- example_gt %>% gt_pseudohaploid(test_n_loci = 3)

# Ploidy is now set to -2
show_ploidy(example_gt)
#> [1] -2

# Individual ploidy now varies between 1 (pseudohaploid) and 2 (diploid)
indiv_ploidy(example_gt)
#> [1] 2 2 1 2 2 1 2
```
